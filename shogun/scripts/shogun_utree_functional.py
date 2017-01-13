"""
Functional repertiore prediction with UTree Species/Strain Classification
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
def shogun_utree_functional(input, output, utree_indx, threads):
    verify_make_dir(output)

    basenames = [os.path.basename(filename)[:-4] for filename in os.listdir(input) if filename.endswith('.fna')]

    for basename in basenames:
        fna_file = os.path.join(input, basename + '.fna')
        tsv_outf = os.path.join(output, basename + 'species.utree.tsv')
        if not os.path.isfile(tsv_outf):
            print(utree_search(utree_indx, fna_file, tsv_outf))
        else:
            print("Found the output file \"%s\". Skipping the alignment phase for this file." % tsv_outf)

    counts = []
    utree_species_outf = os.path.join(output, 'taxon_counts.txt')
    # Indexing for emblalmer
    if not os.path.isfile(utree_species_outf):
        for basename in basenames:
            lcas = []
            utree_tsv = os.path.join(output, basename + 'species.utree.tsv')
            with open(utree_tsv) as inf:
                tsv_parser = csv.reader(inf, delimiter='\t')
                for line in tsv_parser:
                    if line[1]:
                        lcas.append(';'.join(line[1].split('; ')))
            counts.append(Counter(filter(None, lcas)))

    species_df = pd.DataFrame(counts, index=basenames)
    species_df.T.to_csv(utree_species_outf, sep='\t')

if __name__ == '__main__':
    shogun_utree_functional()
