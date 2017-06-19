"""
Copyright 2015-2017 Knights Lab, Regents of the University of Minnesota.

This software is released under the GNU Affero General Public License (AGPL) v3.0 License.
"""

import os
import click

import logging

from shogun.aligners import EmbalmerAligner

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
    'embalmer': EmbalmerAligner
}

@cli.command(help="Run the SHOGUN aligner")
@click.option('-a', '--aligner', type=click.Choice(['bt2', 'embalmer', 'utree']), default='embalmer',
              help='The aligner to use.', show_default=True)
@click.option('-i', '--input', type=click.Path(), required=True, help='The file containing the combined seqs.')
@click.option('-d', '--database', type=click.Path(), required=True, help="The database file.")
@click.option('-o', '--output', type=click.Path(), default=os.getcwd(), help='The output folder directory', show_default=True)
@click.pass_context
def align(ctx, aligner, input, database, output):
    # Set up the logger
    file_handler = logging.FileHandler(os.path.join(output, 'shogun.log'))
    log_formatter = logging.Formatter('%(asctime)s %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')
    file_handler.setFormatter(log_formatter)
    logging.getLogger().addHandler(file_handler)

    aligner = ALIGNERS[aligner](database)

    aligner.align(infile, outfile)





if __name__ == '__main__':
    cli()
