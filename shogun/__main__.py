"""
Copyright 2015-2017 Knights Lab, Regents of the University of Minnesota.

This software is released under the GNU Affero General Public License (AGPL) v3.0 License.
"""

import os
import logging
from datetime import date
from multiprocessing import cpu_count
import click
from yaml import load
import pandas as pd

from shogun import __version__, logger
from shogun.aligners import EmbalmerAligner, UtreeAligner, BowtieAligner
from shogun.coverage import get_coverage_of_microbes
from shogun.function import function_run_and_save, parse_function_db
from shogun.redistribute import redistribute_taxatable, parse_bayes
from shogun.utils import normalize_by_median_depth

ROOT_COMMAND_HELP = """\
SHOGUN command-line interface\n
--------------------------------------
"""


SETTINGS = dict()
TAXA = ['kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species', 'strain']
TAXAMAP = dict(zip(TAXA, range(1, 9)))

@click.group(invoke_without_command=False, help=ROOT_COMMAND_HELP)
@click.option('--log', type=click.Choice(('debug', 'info', 'warning', 'critical')), help="The log level to record.", default="warning")
@click.option('--shell/--no-shell', default=False, help="Use the shell for Python subcommands (not recommended).")
@click.version_option(version=__version__)
@click.pass_context
def cli(ctx, log, shell):
    ctx.obj = {'log': log}

    if log == 'debug':
        logger.setLevel(logging.DEBUG)
    elif log == 'info':
        logger.setLevel(logging.INFO)
    elif log == 'warning':
        logger.setLevel(logging.WARNING)
    else:
        logger.setLevel(logging.CRITICAL)
    ctx.shell = shell


ALIGNERS = {
    'embalmer': EmbalmerAligner,
    'utree': UtreeAligner,
    'bowtie2': BowtieAligner
}

#TODO: Allow for shell=True
#TODO: Turn off post-alignment
#TODO: Fix redistribute bug
# shogun --log debug align --input ./shogun/tests/data/combined_seqs.fna --aligner utree --database /project/flatiron2/analysis_SHOGUN/data/references/rep82 --output ~/scratch_shogun --level off --threads 1
@cli.command(help="Run the SHOGUN aligner")
@click.option('-a', '--aligner', type=click.Choice(['all', 'bowtie2', 'embalmer', 'utree']), default='embalmer',
              help='The aligner to use [Note: default embalmer is capitalist, use embalmer-tax if you want to redistribute].', show_default=True)
@click.option('-i', '--input', type=click.Path(resolve_path=True, exists=True, allow_dash=True), required=True, help='The file containing the combined seqs.')
@click.option('-d', '--database', type=click.Path(resolve_path=True, exists=True), default=os.getcwd(), help="The path to the database folder.")
@click.option('-o', '--output', type=click.Path(resolve_path=True, writable=True), default=os.path.join(os.getcwd(), date.today().strftime('results-%y%m%d')), help='The output folder directory', show_default=True)
@click.option('-l', '--level', type=click.Choice(TAXA + ['all', 'off']), default='strain', help='The level to collapse taxatables and functions to (not required, can specify off).')
@click.option('--function/--no-function', default=True, help='Run functional algorithms. **This will normalize the taxatable by median depth.')
@click.option('--capitalist/--no-capitalist', default=True, help='Run capitalist with embalmer post-align or not.')
@click.option('-t', '--threads', type=click.INT, default=cpu_count(), help="Number of threads to use.")
@click.pass_context
def align(ctx, aligner, input, database, output, level, function, capitalist, threads):
    if not os.path.exists(output):
        os.makedirs(output)

    if not capitalist:
        # Set to not run Embalmer post-align in capitalist mode
        ALIGNERS['embalmer'] = lambda database, threads=threads, shell=ctx.shell: EmbalmerAligner(database, shell=shell, threads=threads, post_align='taxonomy')

    if aligner == 'all':
        redist_outs = []
        redist_levels = []
        for align in ALIGNERS.values():
            aligner_cl = align(database, threads=threads, shell=ctx.shell)
            aligner_cl.align(input, output)
            if level is not 'off':
                redist_out = os.path.join(output, "%s_taxatable.%s.txt" % (aligner_cl._name, level))
                _redist_outs, _redist_levels = _redistribute(database, level, redist_out, aligner_cl.outfile)
                redist_outs.extend(_redist_outs)
                redist_levels.extend(_redist_levels)
    else:
        aligner_cl = ALIGNERS[aligner](database, threads=threads)
        aligner_cl.align(input, output)
        logger.debug(level)
        if level is not 'off':
            redist_out = os.path.join(output, "taxatable.%s.txt" % (level))
            redist_outs, redist_levels = _redistribute(database, level, redist_out, aligner_cl.outfile)

    if function:
        _function(redist_outs, database, output, redist_levels)


@cli.command(help="Run the SHOGUN redistribution algorithm.")
@click.option('-i', '--input', type=click.Path(resolve_path=True, exists=True, allow_dash=True), required=True, help="The taxatable.")
@click.option('-d', '--database', type=click.Path(resolve_path=True, exists=True), required=True, help="The path to the database folder.")
@click.option('-l', '--level', type=click.Choice(TAXA + ['all']), default='strain', help='The level to collapse to.')
@click.option('-o', '--output', type=click.Path(resolve_path=True, writable=True), default=os.path.join(os.getcwd(), date.today().strftime('taxatable-%y%m%d.txt')), help='The output file', show_default=True)
@click.pass_context
def redistribute(ctx, input, database, level, output):
    _redistribute(database, level, output, input)

def _redistribute(database, level, outfile, redist_inf):
    logger.debug("Beginning redistribution for file: %s" % redist_inf)
    data_files = _load_metadata(database)

    shear = os.path.join(database, data_files['general']['shear'])

    shear_df = parse_bayes(shear)

    output_files = []
    output_levels = []

    if level == 'all':
        for l in TAXA:
            df_output = redistribute_taxatable(redist_inf, shear_df, level=TAXAMAP[l])
            tmp_spl = outfile.split('.')
            tmp_path = '.'.join(tmp_spl[:-1] + [l] + [tmp_spl[-1]])
            df_output.to_csv(tmp_path, sep='\t', float_format="%d",na_rep=0, index_label="#OTU ID")
            output_files.append(tmp_path)
            output_levels.append(l)
    else:
        df_output = redistribute_taxatable(redist_inf, shear_df, level=TAXAMAP[level])
        df_output.to_csv(outfile, sep='\t', float_format="%d", na_rep=0, index_label="#OTU ID")
        output_files.append(outfile)
        output_levels.append(level)

    return output_files, output_levels

@cli.command(help="Run the SHOGUN functional algorithm.")
@click.option('-i', '--input', type=click.Path(resolve_path=True, exists=True, allow_dash=True), required=True, help="The taxatable.")
@click.option('-d', '--database', type=click.Path(resolve_path=True, exists=True), required=True, help="The path to the folder containing the function database.")
@click.option('-o', '--output', type=click.Path(resolve_path=True, writable=True), default=os.path.join(os.getcwd(), date.today().strftime('results-%y%m%d')), help='The output file', show_default=True)
@click.option('-l', '--level', type=click.Choice(['genus', 'species', 'strain']), default='strain', help='The level to collapse to.')
@click.pass_context
def function(ctx, input, database, output, level):
    _function([input], database, output, [level])

def _function(inputs, database, output, levels):
    # Check if output exists, if not then make
    if not os.path.exists(output):
        os.makedirs(output)

    # Load the datafiles to locate function db
    data_files = _load_metadata(database)

    # Load the functional db
    logger.info("Loading the functional database and converting.")
    func_db = parse_function_db(data_files, database)

    for input, level in zip(inputs, levels):
        # Verify it is in a reasonable level
        if level in ['genus', 'species', 'strain']:
            logger.info("Starting functional prediction with input file %s at level %s" % (os.path.abspath(input), level))
            function_run_and_save(input, func_db, output, TAXAMAP[level])
        else:
            continue

@cli.command(help="Normalize a taxatable by median depth.")
@click.option('-i', '--input', type=click.Path(resolve_path=True, exists=True, allow_dash=True), required=True, help="The output taxatable.")
@click.option('-o', '--output', type=click.Path(resolve_path=True, writable=True, allow_dash=True), help="The taxatable output normalized by median depth.", default=os.path.join(os.getcwd(), date.today().strftime('taxatable.normalized-%y%m%d.txt')), show_default=True)
def normalize(input, output):
    df = pd.read_csv(input, sep="\t", index_col=0)
    outdf = normalize_by_median_depth(df)
    outdf.to_csv(output, sep='\t', float_format="%d",na_rep=0, index_label="#OTU ID")


@cli.command(help="Show confidence of coverage of microbes.")
@click.option('-i', '--input', type=click.Path(resolve_path=True, exists=True, allow_dash=True), required=True, help="The output BURST capitalist alignment.")
@click.option('-d', '--database', type=click.Path(resolve_path=True, exists=True), required=True, help="The path to the folder containing the function database.")
@click.option('-o', '--output', type=click.Path(resolve_path=True, writable=True), help="The coverage table.", default=os.path.join(os.getcwd(), date.today().strftime('coverage-%y%m%d.txt')), show_default=True)
@click.option('-l', '--level', type=click.Choice(['genus', 'species', 'strain']), default='strain', help='The level to collapse to.')
def coverage(input, database, output, level):
    # This is only the coverage script
    _coverage(input, database, output, TAXAMAP[level])

def _coverage(input, database, output, level):
    data_files = _load_metadata(database)

    shear = os.path.join(database, data_files['general']['shear'])

    shear_df = parse_bayes(shear)
    outdf = get_coverage_of_microbes(input, shear_df, level)
    outdf.to_csv(output, sep='\t', float_format="%.5f",na_rep=0)

def _load_metadata(database):
    metadata_file = os.path.join(database, 'metadata.yaml')
    if os.path.exists(metadata_file):
        with open(metadata_file, 'r') as stream:
            logger.debug(
                "Attempting to load the database metadata file at %s" % (os.path.abspath(metadata_file)))
            data_files = load(stream)
        return data_files
    else:
        logger.critical("Unable to load database at %s" % os.path.abspath(metadata_file))
        raise Exception("Unable to load database at %s" % os.path.abspath(metadata_file))

#TODO: Implement post alignment phase only
def assign():
    pass

if __name__ == '__main__':
    cli()
