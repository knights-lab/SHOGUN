#!/usr/bin/env python
import click
import os
from collections import Counter, defaultdict
import csv
import pandas as pd
import pickle

from ninja_utils.utils import verify_make_dir

from ninja_shogun.wrappers import utree_search


def build_img_map(infile: str):
    gg2img_oid = defaultdict(int)
    df = pd.DataFrame.from_csv(infile)
    for row in df.iterrows():
        gg2img_oid[row] = row
    return gg2img_oid


@click.command()
@click.option('-i', '--input', type=click.Path(), default=os.getcwd(), help='Directory containing the input FASTA files with ".fna" extensions (default=cwd)')
@click.option('-o', '--output', type=click.Path(), default=os.path.join(os.getcwd(), 'shogun_bugbase'), help='Output directory for the results')
@click.option('-u', '--img_database_folder', type=click.Path(), help='Location of the BugBase Database folder.')
def shogun_bugbase(input, output, img_database_folder):
    verify_make_dir(output)
    utree_outf = os.path.join(output, 'taxa_counts.txt')
    # Indexing for emblalmer
    if not os.path.isfile(utree_outf):

        utree_indx = os.path.join(img_database_folder, 'img.genes.ctr')
        with open(os.path.join(img_database_folder, 'img_map.pkl'), 'rb') as inf:
            gg2img_oid = pickle.load(inf)

        basenames = [os.path.basename(filename)[:-4] for filename in os.listdir(input) if filename.endswith('.fna')]

        for basename in basenames:
            fna_file = os.path.join(input, basename + '.fna')
            tsv_outf = os.path.join(output, basename + '.utree.tsv')
            if not os.path.isfile(tsv_outf):
                print(utree_search(utree_indx, fna_file, tsv_outf))
            else:
                print("Found the output file \"%s\". Skipping the alignment phase for this file." % tsv_outf)

        counts = []

        for basename in basenames:
            lcas = []
            utree_tsv = os.path.join(output, basename + '.utree.tsv')
            with open(utree_tsv) as inf:
                tsv_parser = csv.reader(inf, delimiter='\t')
                for line in tsv_parser:
                    if line[1]:
                        taxon = line[1].replace('; ', ';')
                        if taxon in gg2img_oid:
                            lcas.append(gg2img_oid[taxon])
            counts.append(Counter(filter(None, lcas)))

        df = pd.DataFrame(counts, index=basenames).fillna(0).astype(int)
        columns = df.columns
        columns[0] = '#OTU ID'
        df.columns = columns
        df.T.to_csv(utree_outf, sep='\t')
    else:
        print("Found the output file \"%s\". Skipping all steps." % utree_outf)


if __name__ == '__main__':
    shogun_bugbase()
