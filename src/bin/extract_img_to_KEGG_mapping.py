#!/usr/bin/env python
# usage
# python extract_IMG_to_KEGG_mapping.py -i IMGDIR -o outputfile.txt
from __future__ import print_function
import argparse
import re
import csv
from collections import defaultdict
import sys


def make_arg_parser():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('-i', '--IMG_dir', help='The top-level IMG directory', required=True)
    parser.add_argument('-o', '--output', help='If not provided, output to stdout, else write to this file')
    return parser


def read_img_ko_table(f, IMG_ID=None):
    ### if IMG_ID is not None, gene IDs will be returns as IMG_ID_gene_ID
    next(f)

    for line in f:
        words = line.strip().split('\t')
        gene_id = words[0]
        ko_id = words[10]
        if IMG_ID is not None:
            gene_id = IMG_ID + '_' + gene_id
        yield(gene_id, ko_id)


def write_gene_length_counts(outf, lengths):
    columns = ['domain', 'genus', 'species', 'length']

    wr = csv.writer(outf, quoting=csv.QUOTE_ALL)
    wr.writerow(columns)

    map(wr.writerow, [list(species) + [lengths[species]] for species in lengths])
    outf.close()


def main():
    parser = make_arg_parser()
    args = parser.parse_args()

    img_dirs = [f in os.listdir() if os.isdir(f)]

    with open(args., 'rb') as inf:
        database_gen = read_fasta(inf)
        for title, line in database_gen:
            taxon_oid = re.search(find, title).group(0)
            lengths[taxa_map.get(taxon_oid)] += len(line)

    with open(args.output, 'wb') if args.output else sys.stdout as outf:
        write_gene_length_counts(outf, lengths)


if __name__ == '__main__':
    main()
