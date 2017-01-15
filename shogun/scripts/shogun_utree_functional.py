"""
functional repertiore prediction with utree species/strain classification
"""

import click
import os
from collections import counter
import csv
import pandas as pd

from ninja_utils.utils import verify_make_dir

from shogun.wrappers import utree_search


@click.command()
@click.option('-i', '--input', type=click.path(), default=os.getcwd(), help='directory containing the input fasta files with ".fna" extensions (default=cwd)')
@click.option('-o', '--output', type=click.path(), default=os.path.join(os.getcwd(), 'shogun_utree_lca_out'), help='output directory for the results')
@click.option('-u', '--utree_indx', required=true, help='path to the bowtie2 index')
@click.option('-p', '--threads', type=click.int, default=1, help='the number of threads to use (default=1)')
def shogun_utree_functional(input, output, utree_indx, threads):
    verify_make_dir(output)

    basenames = [os.path.basename(filename)[:-4] for filename in os.listdir(input) if filename.endswith('.fna')]

    # Run the species pipeline
    for basename in basenames:
        fna_file = os.path.join(input, basename + '.fna')
        tsv_outf = os.path.join(output, basename + '.species.utree.tsv')
        if not os.path.isfile(tsv_outf):
            print(utree_search(utree_indx, fna_file, tsv_outf))
        else:
            print("found the output file \"%s\". skipping the species alignment phase for this file." % tsv_outf)

    # Run the strain pipeline
    for basename in basenames:
        fna_file = os.path.join(input, basename + '.fna')
        tsv_outf = os.path.join(output, basename + '.strain.utree.tsv')
        if not os.path.isfile(tsv_outf):
            print(utree_search(utree_indx, fna_file, tsv_outf))
        else:
            print("found the output file \"%s\". skipping the strain alignment phase for this file." % tsv_outf)

    for basename in basenames:
        species_tsv = os.path.join(output, basename + '.species.utree.tsv')
        strain_tsv = os.path.join(output, basename + '.strain.utree.tsv')
        with open(species_tsv) as species_inf:
            species_reader = csv.reader(species_inf, delimiter="\t")
            with open(strain_tsv) as strain_inf:
                strain_reader = csv.reader(strain_inf, delimiter="\t")
                for species, strain in zip(species_reader, strain_reader):
                    print(species)
                    print(strain)

if __name__ == '__main__':
    shogun_utree_functional()
