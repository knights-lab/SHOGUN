#!/usr/bin/env python
"""
Copyright 2015-2017 Knights Lab, Regents of the University of Minnesota.

This software is released under the GNU Affero General Public License (AGPL) v3.0 License.
"""

import click
import os
import sys
from collections import Counter
import csv
import re
import pandas as pd
from math import floor
from ninja_utils.utils import verify_make_dir

from shogun.wrappers import utree_search


@click.command()
@click.option('-i', '--input', type=click.Path(), default=os.getcwd(), help='Directory containing the input FASTA files with ".fna" extensions (default=cwd)')
@click.option('-o', '--output', type=click.Path(), default=os.path.join(os.getcwd(), 'shogun_utree_lca_out'), help='Output directory for the results')
@click.option('-u', '--utree_indx', required=True, help='Path to the bowtie2 index')
@click.option('-p', '--threads', type=click.INT, default=1, help='The number of threads to use (default=1)')
@click.option('-c', '--confidence', type=click.FLOAT, default=.8, help='Required confidence threshold for kmer matching calculated as 1 - <support for second best> / <support for best> (0 to 1.0) (default=0.8)')
@click.option('-c', '--support', type=click.INT, default=5, help='Minimum number of sliding-window kmers spaced at least 4 bases apart throughout the query matched the most-matched species? (1 to ~25, depends on query lengths) (default=5)')
@click.option('-m', '--mincount', type=click.INT, default=20, help='Minimum count (number of reads matching) per species (default=20)')

def shogun_utree_lca(input, output, utree_indx, threads, confidence, support, mincount):
    verify_make_dir(output)

    basenames = [os.path.basename(filename)[:-4] for filename in os.listdir(input) if filename.endswith('.fna')]

    for basename in basenames:
        fna_file = os.path.join(input, basename + '.fna')
        tsv_outf = os.path.join(output, basename + '.utree.tsv')
        if not os.path.isfile(tsv_outf):
            print(utree_search(utree_indx, fna_file, tsv_outf))
        else:
            print("Found the output file \"%s\". Skipping the alignment phase for this file." % tsv_outf)

    counts = []
    utree_outf = os.path.join(output, 'taxon_counts.txt')
    # Tabulating
    print("Tabulating and filtering hits...")
    
    # print a row of "-" for every 10 samples
    if len(basenames) >= 100:
        for i in range(floor(len(basenames)/10)):
            sys.stdout.write('-')
        sys.stdout.write('\n')
        sys.stdout.flush()
    if not os.path.isfile(utree_outf):
        n_fail_confidence_only = 0
        n_fail_support_only = 0
        n_fail_both = 0
        n = 0
        n_pass = 0
        for i,basename in enumerate(basenames):
            if len(basenames) >= 100:
                if (i+1) % 10 == 0:
                    sys.stdout.write('.')
                    sys.stdout.flush()
            lcas = [] # list of tuples [taxonomy, confidence, support]
            utree_tsv = os.path.join(output, basename + '.utree.tsv')
            with open(utree_tsv) as inf:
                tsv_parser = csv.reader(inf, delimiter='\t')
                for line in tsv_parser:
                    if line[1]:
                        taxonomy  = line[1]
                        is_confident = float(line[2]) >= confidence
                        is_supported = int(line[3]) >= support
                        n += 1
                        if not is_confident and not is_supported:
                            n_fail_both += 1
                        elif not is_confident:
                            n_fail_confidence_only += 1
                        elif not is_supported:
                            n_fail_support_only += 1
                        else:
                            n_pass += 1
                            lcas.append(taxonomy)
            counts.append(Counter(lcas))
        print('%d total assignments\n%d failed confidence only\n%d failed support_only\n%d failed both\n%d remaining' %(n,n_fail_confidence_only,n_fail_support_only,n_fail_both,n_pass))
    sys.stdout.write('\n')
    sys.stdout.flush()

    df = pd.DataFrame(counts, index=basenames)
    # filter by mincount
    df[df < mincount] = 0
    # drop spaces in column
    df.columns = [colname.replace('; ',';') for colname in df.columns]
    # drop trailing t__ in taxonomy
    df.columns = [re.sub(';t__$','',colname) for colname in df.columns]
    df.T.to_csv(os.path.join(output, 'taxon_counts.csv'),
                index_label='Taxon',na_rep='0',sep='\t')

if __name__ == '__main__':
    shogun_utree_lca()
