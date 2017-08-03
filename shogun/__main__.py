"""
Copyright 2015-2017 Knights Lab, Regents of the University of Minnesota.

This software is released under the GNU Affero General Public License (AGPL) v3.0 License.
"""

import logging
import os
from datetime import date
from multiprocessing import cpu_count

import click
from yaml import load

from shogun import __version__
from shogun.aligners import EmbalmerAligner, UtreeAligner, BowtieAligner
from shogun.function import function_run_and_save, parse_function_db
from shogun.redistribute import redistribute_taxatable, parse_bayes

ROOT_COMMAND_HELP = """\
SHOGUN command-line interface\n
--------------------------------------
"""


SETTINGS = dict()
TAXA = ['kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species', 'strain']
TAXAMAP = dict(zip(TAXA, range(1, 9)))

@click.group(invoke_without_command=False, help=ROOT_COMMAND_HELP)
@click.option('--debug/--no-debug', default=False)
@click.option('--verbose/--no-verbose', default=True)
@click.version_option(version=__version__)
@click.pass_context
def cli(ctx, debug, verbose):
    ctx.obj = {'DEBUG': debug, 'VERBOSE': verbose}
    # Setup the logger
    log_formatter = logging.Formatter('%(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')
    root_logger = logging.getLogger()

    if ctx.obj['DEBUG']:
        root_logger.setLevel(logging.DEBUG)
    elif verbose:
        root_logger.setLevel(logging.INFO)
    else:
        root_logger.setLevel(logging.WARNING)

    console_handler = logging.StreamHandler()
    console_handler.setFormatter(log_formatter)
    root_logger.addHandler(console_handler)
    global root_logger


ALIGNERS = {
    'embalmer': EmbalmerAligner,
    'utree': UtreeAligner,
    'bowtie2': BowtieAligner
}

@cli.command(help="Run the SHOGUN aligner")
@click.option('-a', '--aligner', type=click.Choice(['all', 'bowtie2', 'embalmer', 'utree']), default='embalmer',
              help='The aligner to use.', show_default=True)
@click.option('-i', '--input', type=click.Path(), required=True, help='The file containing the combined seqs.')
@click.option('-d', '--database', type=click.Path(), default=os.getcwd(), help="The path to the database folder.")
@click.option('-o', '--output', type=click.Path(), default=os.path.join(os.getcwd(), date.today().strftime('results-%y%m%d')), help='The output folder directory', show_default=True)
@click.option('-l', '--level', type=click.Choice(TAXA + ['all', 'off']), default='strain', help='The level to collapse taxatables and functions to (not required, can specify off).')
@click.option('--function/--no-function', default=True, help='Run functional algorithms.')
@click.option('-t', '--threads', type=click.INT, default=cpu_count(), help="Number of threads to use.")
@click.pass_context
def align(ctx, aligner, input, database, output, level, function, threads):
    if not os.path.exists(output):
        os.makedirs(output)

    # Set up the logger
    file_handler = logging.FileHandler(os.path.join(output, 'shogun.log'))
    log_formatter = logging.Formatter('%(asctime)s %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')
    file_handler.setFormatter(log_formatter)
    logging.getLogger().addHandler(file_handler)

    if aligner == 'all':
        redist_outs = []
        redist_levels = []
        for align in ALIGNERS.values():
            aligner_cl = align(database, threads=threads)
            aligner_cl.align(input, output)
            if level is not 'off':
                redist_out = os.path.join(output, "%s_taxatable.%s.txt" % (aligner_cl._name, level))
                _redist_outs, _redist_levels = _redistribute(database, level, redist_out, aligner_cl.outfile)
                redist_outs.extend(_redist_outs)
                redist_levels.extend(_redist_levels)
    else:
        aligner_cl = ALIGNERS[aligner](database, threads=threads)
        aligner_cl.align(input, output)
        if level is not 'off':
            redist_out = os.path.join(output, "taxatable.%s.txt" % (level))
            redist_outs, redist_levels = _redistribute(database, level, redist_out, aligner_cl.outfile)

    if function:
        _function(redist_outs, database, output, redist_levels)


@cli.command(help="Run the SHOGUN redistribution algorithm.")
@click.option('-i', '--input', type=click.Path(), required=True, help="The taxatable.")
@click.option('-d', '--database', type=click.Path(), required=True, help="The path to the database folder.")
@click.option('-l', '--level', type=click.Choice(TAXA + ['all']), default='strain', help='The level to collapse to.')
@click.option('-o', '--output', type=click.Path(), default=os.path.join(os.getcwd(), date.today().strftime('taxatable-%y%m%d.txt')), help='The output file', show_default=True)
@click.pass_context
def redistribute(ctx, input, database, level, output):
    _redistribute(database, level, output, input)

def _redistribute(database, level, outfile, redist_inf):
    with open(os.path.join(database, 'metadata.yaml'), 'r') as stream:
        data_files = load(stream)

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
@click.option('-i', '--input', type=click.Path(), required=True, help="The taxatable.")
@click.option('-d', '--database', type=click.Path(), required=True, help="The path to the folder containing the function database.")
@click.option('-o', '--output', type=click.Path(), default=os.path.join(os.getcwd(), date.today().strftime('results-%y%m%d')), help='The output file', show_default=True)
@click.option('-l', '--level', type=click.Choice(['genus', 'species', 'strain']), default='strain', help='The level to collapse to.')
@click.pass_context
def function(ctx, input, database, output, level):
    _function([input], database, output, [level])

def _function(inputs, database, output, levels):
    # Check if output exists, if not then make
    if not os.path.exists(output):
        os.makedirs(output)

    with open(os.path.join(database, 'metadata.yaml'), 'r') as stream:
        data_files = load(stream)

    func_db = parse_function_db(data_files, database)

    for input, level in zip(inputs, levels):
        # Verify it is in a reasonable level
        if level in ['genus', 'species', 'strain']:
            root_logger.info("Starting functional prediction with input file %s at level %s" % (input, level))
            function_run_and_save(input, func_db, output, TAXAMAP[level])
        else:
            continue

if __name__ == '__main__':
    cli()
