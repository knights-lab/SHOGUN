#!/usr/bin/env python
"""
Copyright 2015-2017 Knights Lab, Regents of the University of Minnesota.

This software is released under the GNU Affero General Public License (AGPL) v3.0 License.
"""

import click
from collections import defaultdict
import os

from ninja_utils.parsers import FASTA

@click.command()
@click.option('-i', '--input', type=click.File(), default='-', help='The input annotated FASTA file')
@click.option('-m', '--map', type=click.File(), default='-', help='The input FASTA mapping file for UTree redistribute')
@click.option('-o', '--output', type=click.File('w'), default=os.path.join(os.getcwd(), 'annotated'), help='The directory to output the formatted DB and BT2 db (default=annotated)')
def extract_genome_lengths(input, map, output):
    d = defaultdict(int)
    inf_fasta = FASTA(input)
    for header, seq in inf_fasta.read():
        map_line = next(map).rstrip().split('\t')[1].replace('; ', ';')
        d[map_line] += len(seq)
    for key, value in d.items():
        output.write('%s\t%s\n' % (key, value))

if __name__ == '__main__':
    extract_genome_lengths()
