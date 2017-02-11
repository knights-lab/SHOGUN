#!/usr/bin/env python
"""
Copyright 2015-2017 Knights Lab, Regents of the University of Minnesota.

This software is released under the GNU Affero General Public License (AGPL) v3.0 License.
"""

import click
import os
from collections import Counter
import csv
import pandas as pd

from ninja_utils.utils import verify_make_dir

from shogun.wrappers import utree_search


@click.command()
@click.option('-i', '--input', type=click.Path(), default=os.getcwd(), help='Directory containing the input FASTA files with ".fna" extensions (default=cwd)')
@click.option('-o', '--output', type=click.Path(), default=os.path.join(os.getcwd(), 'shogun_utree_lca_out'), help='Output directory for the results')
@click.option('-u', '--utree_indx', required=True, help='Path to the bowtie2 index')
@click.option('-p', '--threads', type=click.INT, default=1, help='The number of threads to use (default=1)')
@click.option('-c', '--confidence', type=click.FLOAT, default=.8, help='Required confidence threshold for kmer matching calculated as 1 - <support for second best> / <support for best> (0 to 1.0) (default=0.8)')
@click.option('-c', '--support', type=click.INT, default=5, help='Minimum number of sliding-window kmers spaced at least 4 bases apart throughout the query matched the most-matched species? (1 to ~25, depends on query lengths) (default=5)')
@click.option('-m', '--mincover', type=click.INT, default=20, help='Minimum coverage (number of reads matching) per species (default=20)')

def shogun_utree_lca(input, output, utree_indx, threads, confidence, support, mincover):
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
    # Indexing for emblalmer
    if not os.path.isfile(utree_outf):
        n_fail_confidence_only = 0
        n_fail_support_only = 0
        n_fail_both = 0
        n = 0
        n_remain = 0
        for basename in basenames:
            lcas = [] # list of tuples [taxonomy, confidence, support]
            utree_tsv = os.path.join(output, basename + '.utree.tsv')
            with open(utree_tsv) as inf:
                tsv_parser = csv.reader(inf, delimiter='\t')
                for line in tsv_parser:
                    if line[1]:
                        taxonomy  = ';'.join(line[1].split('; '))
                        n += 1
                        if float(line[2]) < confidence and int(line[3]) < support:
                            n_fail_both += 1
                        elif float(line[2]) < confidence:
                            n_fail_confidence_only += 1
                        elif int(line[3]) < support:
                            n_fail_support_only += 1
                        else:
                            n_remain += 1
                            lcas.append(taxonomy)
            c=Counter(filter(None, lcas))
            c = Counter(el for el in c.elements() if c[el] > mincover)
            counts.append(c)
        print('%d total assignments\n%d failed confidence only\n%d failed support_only\n%d failed both\n%d remaining' %(n,n_fail_confidence_only,n_fail_support_only,n_fail_both,n_remain))

    df = pd.DataFrame(counts, index=basenames)
    df.T.to_csv(os.path.join(output, 'taxon_counts.csv'))

if __name__ == '__main__':
    shogun_utree_lca()
