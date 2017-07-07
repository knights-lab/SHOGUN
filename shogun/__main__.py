"""
Copyright 2015-2017 Knights Lab, Regents of the University of Minnesota.

This software is released under the GNU Affero General Public License (AGPL) v3.0 License.
"""

import os
import click
from datetime import date
import logging
from yaml import load
import glob.glob as glob

from shogun.aligners import EmbalmerAligner, UtreeAligner, BowtieAligner
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

def _check_function_db(metadata: dict, database: str) -> bool:
    if not 'function' in metadata:
        return False
    else:
        file_set = set(glob(os.path.join(database, metadata['function'] + '*')))
        suffices = ['module-annotations.txt', 'species2ko.txt', 'strain2ko.txt']
        files = ["%s-%s" % (os.path.join(database, metadata['function']), suffix) for suffix in suffices]
        for file in files:
            if file not in file_set:
                return False
        return dict(zip(('modules', 'species', 'strain'), files))

@cli.command(help="Run the SHOGUN redistrbution algorithm.")
@click.option('-i', '--input', type=click.Path(), required=True, help="The the taxatable.")
@click.option('-c', '--function', type=click.Path(), required=True, help="The path to the folder containing the functional database.")
@click.option('-o', '--output', type=click.Path(), default=os.path.join(os.getcwd(), date.today().strftime('taxatable-%y%m%d.txt')), help='The output file', show_default=True)
@click.pass_context
def function(ctx, input, function, output):
    pass

def _function(input, function, output):
    pass

if __name__ == '__main__':
    cli()
