"""
Copyright 2015-2017 Knights Lab, Regents of the University of Minnesota.

This software is released under the GNU Affero General Public License (AGPL) v3.0 License.
"""

import os
import click

ROOT_COMMAND_HELP = """\
SHOGUN command-line interface\n
--------------------------------------
"""


@click.group(invoke_without_command=False, help=ROOT_COMMAND_HELP)
@click.option('--debug/--no-debug', default=False)
@click.pass_context
def cli(ctx, debug):
    ctx.obj['DEBUG'] = debug
    # TODO: Setup the logger
    # ctx.obj['LOGGER'] = Logger(logfp=SETTINGS.settings['log'], log_persist=SETTINGS.settings['log_persists'])


@cli.command(help="Run the SHOGUN aligner")
@click.option('-a', '--aligner', type=click.Choice(['bt2', 'embalmer', 'utree']), default='embalmer',
              help='The aligner to use.', show_default=True)
@click.pass_context
def align(ctx):
    click.echo('Debug is %s' % (ctx.obj['DEBUG'] and 'on' or 'off'))


if __name__ == '__main__':
    cli(obj={})
