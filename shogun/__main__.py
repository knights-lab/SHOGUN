"""
Copyright 2015-2017 Knights Lab, Regents of the University of Minnesota.

This software is released under the GNU Affero General Public License (AGPL) v3.0 License.
"""

import os
import click
from datetime import date
import logging

from shogun.aligners import EmbalmerAligner, UtreeAligner, BowtieAligner

ROOT_COMMAND_HELP = """\
SHOGUN command-line interface\n
--------------------------------------
"""

SETTINGS = dict()


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
@click.pass_context
def align(ctx, aligner, input, database, output):
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
    else:
        aligner_cl = ALIGNERS[aligner](database)
        aligner_cl.align(input, output)



if __name__ == '__main__':
    cli()
