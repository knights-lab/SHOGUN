"""
Copyright 2015-2017 Knights Lab, Regents of the University of Minnesota.

This software is released under the GNU Affero General Public License (AGPL) v3.0 License.
"""

import os
import click
from datetime import date
import logging
from yaml import load

from shogun.aligners import EmbalmerAligner, UtreeAligner, BowtieAligner
from shogun.taxonomy import pie_chart_taxatable

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
@click.option('-l', '--level', type=click.Choice(TAXA + ['off']), default='strain', help='The level to collapse too (not required, can specify off).')
@click.pass_context
def align(ctx, aligner, input, database, output, level):
    if not os.path.exists(output):
        os.makedirs(output)

    # Set up the logger
    file_handler = logging.FileHandler(os.path.join(output, 'shogun.log'))
    log_formatter = logging.Formatter('%(asctime)s %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')
    file_handler.setFormatter(log_formatter)
    logging.getLogger().addHandler(file_handler)

    if aligner == 'all':
        for align in ALIGNERS.values():
            aligner_cl = align(database)
            aligner_cl.align(input, output)
            if level is not 'off':
                with open(os.path.join(database, 'metadata.yaml'), 'r') as stream:
                    data_files = load(stream)
                shear = data_files['shear']
                redist_inf = os.path.join(aligner_cl.outfile, "%s_taxatable.%s.txt" % (aligner_cl._name, level))
                _redistribute(input, shear, level, redist_inf)
    else:
        aligner_cl = ALIGNERS[aligner](database)
        aligner_cl.align(input, output)
        if level is not 'off':
            with open(os.path.join(database, 'metadata.yaml'), 'r') as stream:
                data_files = load(stream)
            shear = data_files['shear']
            redist_inf = os.path.join(aligner_cl.outfile, "taxatable.%s.txt" % (level))
            _redistribute(input, shear, level, redist_inf)


@cli.command(help="Run the SHOGUN aligner")
@click.option('-i', '--input', type=click.Path(), required=True, help="The the taxatable.")
@click.option('-s', '--shear', type=click.Path(), required=True, help="The path to the sheared results.")
@click.option('-l', '--level', type=click.Choice(TAXA), default='strain', help='The level to collapse too.')
@click.option('-o', '--output', type=click.File(), default=os.path.join(os.getcwd(), date.today().strftime('taxatable-%y%m%d.txt')), help='The output folder directory', show_default=True)
@click.pass_context
def redistribute(input, shear, level, output):
    _redistribute(shear, level, output, input)


def _redistribute(shear, level, outfile, redist_inf):
    df_output = pie_chart_taxatable(shear, redist_inf, level=TAXAMAP[level])
    df_output.to_csv(outfile, sep='\t', float_format="%d",na_rep=0, index_label="#OTU ID")

if __name__ == '__main__':
    cli()
