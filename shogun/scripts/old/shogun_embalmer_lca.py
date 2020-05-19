#!/usr/bin/env python
"""
Copyright 2015-2020 Knights Lab, Regents of the University of Minnesota.

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

from shogun.wrappers import embalmer_search


@click.command()
@click.option('-i', '--input_dir', type=click.Path(), default=None, help='Input folder containing multiple fasta files, one per sample, with ".fna" or ".fasta" extensions [required]')
@click.option('-o', '--output_dir', type=click.Path(), default=None, help='Output directory for the results (default input directory)')
@click.option('-r', '--embalmer_db', required=True, help='Path to the embalmer database basename. There must be 3 files with this basename and .edb, .acc, .tax extensions.')
@click.option('-p', '--threads', type=click.INT, default=1, help='The number of threads to use (default=1)')
@click.option('-d', '--pct_id', type=click.FLOAT, default=.97, help='The percent ID for alignments (default=.97)')
@click.option('-m', '--mincount', type=click.INT, default=2, help='Minimum count (number of reads matching) per taxon (default=2)')
@click.option('-t', '--taxa_ncbi', default=True, is_flag=True, help='Pass --taxa_ncbi to embalmer')
def shogun_embalmer_lca(input_dir, output_dir, embalmer_db, threads, pct_id, mincount, taxa_ncbi):
    if output_dir is None:
        output_dir = input_dir
    verify_make_dir(output_dir)

    inputfiles = [filename for filename in os.listdir(input_dir) if filename.endswith('.fna') or filename.endswith('fasta')]
    basenames = [os.path.splitext(filename)[0] for filename in inputfiles]
    outputfps = []

    for i,filename in enumerate(inputfiles):
        input_fp = os.path.join(input_dir, filename)
        tsv_outf = os.path.join(output_dir, basenames[i] + '.embalmer.tsv')
        outputfps.append(tsv_outf)
        if not os.path.isfile(tsv_outf):
            print("Did not file the output file \"%s\". Running the alignment phase for this file." % tsv_outf)
            print(embalmer_search(input_fp, tsv_outf, embalmer_db+".edb", embalmer_db+".tax", embalmer_db+".acc", threads, pct_id, taxa_ncbi))
        else:
            print("Found the output file \"%s\". Skipping the alignment phase for this file." % tsv_outf)

    counts = []

    # Tabulating
    print("Tabulating and filtering hits...")

    # print a row of "-" for every 10 samples
    if len(inputfiles) >= 100:
        for i in range(floor(len(basenames)/10)):
            sys.stdout.write('-')
        sys.stdout.write('\n')
        sys.stdout.flush()

    taxon_outf = os.path.join(output_dir, 'taxon_counts.tsv')
    if os.path.isfile(taxon_outf):
        print("Skipping tabulation step, output file %s already exists." %(taxon_outf))
    else:
        for outputfp in outputfps:
            with open(outputfp) as output_file:
                tsv_parser = csv.reader(output_file, delimiter='\t')
                taxon_counts = Counter()
                for line in tsv_parser:
                    taxon = line[12]
                    # drop trailing t__ in redistribute
                    taxon = re.sub('; t__$','',taxon)
                    taxon = re.sub('; t__None$','',taxon)
                    taxon_counts[taxon] += 1
            counts.append(taxon_counts)
        df = pd.DataFrame(counts, index=basenames)
        # filter by mincount
        df[df < mincount] = 0
        # drop spaces in column
        df.columns = [colname.replace('; ',';') for colname in df.columns]
        # drop columns that sum to zero
        df = df.loc[:,(df.sum(axis=0) != 0)]
        df.T.to_csv(taxon_outf,
                index_label='Taxon',na_rep='0',sep='\t')

        get_rank_specific_taxonomy_tables(df,output_dir)


def get_rank_specific_taxonomy_tables(df, output_dir, extrapolate=True):
    taxa = df.columns

    tables = [] # will be a list of dataframes

    # aggregate "up" -- add descendant counts to each taxon
    for level in range(8):
        taxa_level = [';'.join(taxon.split(';')[:(level+1)]) for taxon in taxa if len(taxon.split(';')) >= level+1]
        unique_taxa = sorted(set(taxa_level))
        # new data frame for this level full of zeros
        df_i = pd.DataFrame(0, index=df.index, columns=unique_taxa)
        # fill one taxon at a time
        for taxon in unique_taxa:
            df_i[taxon] = df.loc[:,[column.startswith(taxon) for column in df.columns]].sum(axis=1)
        tables.append(df_i)

    # interpolate "down" -- distribute parent counts to descendants
    tables_norm = [tables[0].copy()]
    for level in range(1,8):
        df_i = tables[level].copy()
        unique_taxa = tables[level-1].columns
        for taxon in unique_taxa:
            children_ix = [column.startswith(taxon) for column in tables[level]]
            child_sum = df_i.iloc[:,children_ix].sum(axis=1) # row sums
            weights = df_i.iloc[:,children_ix].apply(lambda col: col / child_sum, axis=0)
            new_child_counts = weights.apply(lambda col: col * tables_norm[level-1][taxon], axis=0)
            df_i.iloc[:,children_ix] = new_child_counts
        tables_norm.append(df_i)

    # write files
    for level in range(len(tables)):
        outf = os.path.join(output_dir, 'taxon_counts_L%d.tsv' %(level+1))
        tables[level].T.to_csv(outf,
                index_label='Taxon',na_rep='0', sep='\t')
        outf = os.path.join(output_dir, 'taxon_counts_interpolate_L%d.tsv' %(level+1))
        tables_norm[level].T.to_csv(outf,
                index_label='Taxon',na_rep='0', sep='\t')


    # if extrapolate, then project counts down



if __name__ == '__main__':
    shogun_embalmer_lca()

