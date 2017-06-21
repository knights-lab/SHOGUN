#!/usr/bin/env python
"""
Copyright 2015-2017 Knights Lab, Regents of the University of Minnesota.

This software is released under the GNU Affero General Public License (AGPL) v3.0 License.
"""

import click

from shogun.wrappers import bowtie2_align


@click.command()
@click.argument('infile', type=click.Path(exists=True))
@click.argument('outfile', type=click.Path(exists=False))
@click.argument('database', type=click.Path(exists=True))
def shogun_bt2_align(infile, outfile, database):
    print(bowtie2_align(infile, outfile, database))

if __name__ == '__main__':
    shogun_bt2_align()
