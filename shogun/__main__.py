"""
Copyright 2015-2017 Knights Lab, Regents of the University of Minnesota.

This software is released under the GNU Affero General Public License (AGPL) v3.0 License.
"""

import os
import click
from datetime import date
import logging
from yaml import load
import glob
import csv
from collections import Counter, defaultdict

from scipy.sparse import csr_matrix
import numpy as np
import pandas as pd

from shogun.aligners import EmbalmerAligner, UtreeAligner, BowtieAligner
from shogun.utils import save_csr_matrix, load_csr_matrix
from shogun.taxonomy import pie_chart_taxatable, parse_bayes
from multiprocessing import cpu_count

ROOT_COMMAND_HELP = """\
SHOGUN command-line interface\n
--------------------------------------
"""

SETTINGS = dict()
TAXA = ['kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species', 'strain']
TAXAMAP = dict(zip(TAXA, range(1, 9)))

@click.group(invoke_without_command=False, help=ROOT_COMMAND_HELP)
@click.option('--debug/--no-debug', default=False)
@click.pass_context
def cli(ctx, debug):
    ctx.obj = {'DEBUG': debug}
    # Setup the logger
    log_formatter = logging.Formatter('%(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')
    root_logger = logging.getLogger()

    if ctx.obj['DEBUG']:
        root_logger.setLevel(logging.DEBUG)
    else:
        root_logger.setLevel(logging.INFO)

    console_handler = logging.StreamHandler()
    console_handler.setFormatter(log_formatter)
    root_logger.addHandler(console_handler)


ALIGNERS = {
    'embalmer': EmbalmerAligner,
    'utree': UtreeAligner,
    'bowtie2': BowtieAligner
}

@cli.command(help="Run the SHOGUN aligner")
@click.option('-a', '--aligner', type=click.Choice(['all', 'bowtie2', 'embalmer', 'utree']), default='embalmer',
              help='The aligner to use.', show_default=True)
@click.option('-i', '--input', type=click.Path(), required=True, help='The file containing the combined seqs.')
@click.option('-d', '--database', type=click.Path(), default=os.getcwd(), help="The database file.")
@click.option('-o', '--output', type=click.Path(), default=os.path.join(os.getcwd(), date.today().strftime('results-%y%m%d')), help='The output folder directory', show_default=True)
@click.option('-l', '--level', type=click.Choice(TAXA + ['all', 'off']), default='strain', help='The level to collapse too (not required, can specify off).')
@click.option('-t', '--threads', type=click.INT, default=cpu_count(), help="Number of threads to use.")
@click.pass_context
def align(ctx, aligner, input, database, output, level, threads):
    if not os.path.exists(output):
        os.makedirs(output)

    # Set up the logger
    file_handler = logging.FileHandler(os.path.join(output, 'shogun.log'))
    log_formatter = logging.Formatter('%(asctime)s %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')
    file_handler.setFormatter(log_formatter)
    logging.getLogger().addHandler(file_handler)

    if aligner == 'all':
        for align in ALIGNERS.values():
            aligner_cl = align(database, threads=threads)
            aligner_cl.align(input, output)
            if level is not 'off':
                with open(os.path.join(database, 'metadata.yaml'), 'r') as stream:
                    data_files = load(stream)
                shear = os.path.join(database, data_files['general']['shear'])
                redist_out = os.path.join(output, "%s_taxatable.%s.txt" % (aligner_cl._name, level))
                _redistribute(shear, level, redist_out, aligner_cl.outfile)
    else:
        aligner_cl = ALIGNERS[aligner](database, threads=threads)
        aligner_cl.align(input, output)
        if level is not 'off':
            with open(os.path.join(database, 'metadata.yaml'), 'r') as stream:

                data_files = load(stream)
            shear = os.path.join(database, data_files['general']['shear'])
            redist_out = os.path.join(output, "taxatable.%s.txt" % (level))
            _redistribute(shear, level, redist_out, aligner_cl.outfile)


@cli.command(help="Run the SHOGUN redistribution algorithm.")
@click.option('-i', '--input', type=click.Path(), required=True, help="The taxatable.")
@click.option('-s', '--shear', type=click.Path(), required=True, help="The path to the sheared results.")
@click.option('-l', '--level', type=click.Choice(TAXA + ['all']), default='strain', help='The level to collapse too.')
@click.option('-o', '--output', type=click.Path(), default=os.path.join(os.getcwd(), date.today().strftime('taxatable-%y%m%d.txt')), help='The output file', show_default=True)
@click.pass_context
def redistribute(ctx, input, shear, level, output):
    _redistribute(shear, level, output, input)


def _redistribute(shear, level, outfile, redist_inf):
    shear_df = parse_bayes(shear)
    print(shear, level, outfile, redist_inf)
    if level == 'all':
        for l in TAXA:
            df_output = pie_chart_taxatable(redist_inf, shear_df, level=TAXAMAP[l])
            tmp_spl = outfile.split('.')
            tmp_path = '.'.join(tmp_spl[:-1] + [l] + [tmp_spl[-1]])
            df_output.to_csv(tmp_path, sep='\t', float_format="%d",na_rep=0, index_label="#OTU ID")
    else:
        df_output = pie_chart_taxatable(redist_inf, shear_df, level=TAXAMAP[level])
        df_output.to_csv(outfile, sep='\t', float_format="%d",na_rep=0, index_label="#OTU ID")

def _parse_function_db(metadata: dict, database: str) -> dict:
    if not 'function' in metadata:
        return {}
    else:
        file_set = set(glob.glob(os.path.join(database, metadata['function'] + '*')))
        suffices = ['module-annotations.txt', 'strain2ko.txt']
        files = ["%s-%s" % (os.path.join(database, metadata['function']), suffix) for suffix in suffices]
        for file in files:
            if file not in file_set:
                return {}

        modules_df = _parse_modules(files[0])
        #TODO: Implement the save csr, this works but requires write permissions to db folder
        npz = os.path.join(database, metadata['function']) + "-strain2ko.npz"
        # if not npz in file_set:
        #    row_names, column_names, csr = _parse_kegg_table(files[1])
        #    save_csr_matrix(npz, csr, row_names, column_names)
        # else:
        #   row_names, column_names, csr = load_csr_matrix(npz)
        row_names, column_names, csr = _parse_kegg_table(files[1])
        return dict(zip(('modules_file', 'strain_file', 'strain_names', 'kegg_ids', 'csr', 'modules'), files + [row_names, column_names, csr, modules_df]))

def _parse_modules(infile):
    modules_keggs = defaultdict(Counter)
    with open(infile) as inf:
        csv_inf = csv.reader(inf, delimiter="\t")
        for row in csv_inf:
            modules_keggs[row[0]].update([row[-1][:7]])
    return pd.DataFrame(modules_keggs).fillna(0.0).astype(int)


def _parse_kegg_table(infile):
    indptr = [0]
    indices = []
    kegg_ids = {}
    row_names = {}
    data = []
    with open(infile) as inf:
        csv_inf = csv.reader(inf, delimiter="\t")
        for line in csv_inf:
            row_names.setdefault(line[0], len(row_names))
            counts = Counter(line[1:])
            for key, value in counts.items():
                indices.append(kegg_ids.setdefault(key, len(kegg_ids)))
                data.append(value)
            indptr.append(len(indices))
    return row_names, kegg_ids, csr_matrix((data, indices, indptr), dtype=np.int8)

@cli.command(help="Run the SHOGUN redistribution algorithm.")
@click.option('-i', '--input', type=click.Path(), required=True, help="The the taxatable.")
@click.option('-d', '--database', type=click.Path(), required=True, help="The path to the folder containing the function database.")
@click.option('-o', '--output', type=click.Path(), default=os.path.join(os.getcwd(), date.today().strftime('taxatable-%y%m%d.txt')), help='The output file', show_default=True)
@click.option('-l', '--level', type=click.Choice(['species', 'strain']), default='strain', help='The level to collapse too.')
@click.pass_context
def function(ctx, input, database, output, level):
    _prep_and_do_functions(input, database, output, TAXAMAP[level])

def _prep_and_do_functions(input, database, output, level):
    with open(os.path.join(database, 'metadata.yaml'), 'r') as stream:
        data_files = load(stream)

    db = _parse_function_db(data_files, database)

    kegg_modules_df = db['modules']
    strain_names = db['strain_names']
    kegg_ids = db['kegg_ids']
    kegg_table_csr = db['csr']

    taxatable_df = pd.read_csv(input, sep="\t", index_col=0)
    taxatable_df = taxatable_df[[type(_) == str for _ in taxatable_df.index]]

    taxatable_df['summary'] = [';'.join(_.split(';')[:level]) for _ in taxatable_df.index]
    taxatable_df = taxatable_df.groupby('summary').sum()

    function_df, kegg_modules = _do_function(taxatable_df, strain_names, kegg_ids, kegg_table_csr, kegg_modules_df)



def _do_function(taxatable_df, row_names, column_names, kegg_table_csr, kegg_modules_df):
    _, num_kegg_ids = kegg_table_csr.shape
    _, num_samples = taxatable_df.shape

    kegg_table = np.zeros((num_samples, num_kegg_ids), dtype=np.int)
    for i, row in taxatable_df.iterrows():
        if row.name in row_names:
            idx = row_names[row.name]
            kegg_table += np.outer(row, kegg_table_csr.getrow(idx).todense())

    out_kegg_table_df = pd.DataFrame(kegg_table, index=taxatable_df.columns, columns=sorted(column_names, key=column_names.get), dtype=np.int)
    # kegg modules df
    out_kegg_modules_df = kegg_modules_df.dot(out_kegg_table_df.loc[:,kegg_modules_df.columns].T.fillna(0))
    return out_kegg_table_df, out_kegg_modules_df

if __name__ == '__main__':
    cli()
