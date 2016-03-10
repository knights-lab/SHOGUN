#!/usr/bin/env python
# usage
# python extract_IMG_to_KEGG_mapping.py -i IMGDIR -o outputfile.txt
from __future__ import print_function
import argparse
import csv
import sys
import os

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
        ko_id = words[9]
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

    # list all IMG strain directories
    sys.stderr.write('Listing all IMG directories...\n')
    sys.stderr.flush()
    img_dirs = [f for f in os.listdir(args.IMG_dir) if os.path.isdir(os.path.join(args.IMG_dir, f))]

    sys.stderr.write('Parsing %d directories...\n' %(len(img_dirs)))
    sys.stderr.flush()

    with open(args.output, 'wb') if args.output else sys.stdout as outf:
        for i, imgID in enumerate(img_dirs):
            if (i + 1) % 100 == 0:
                sys.stderr.write('Parsing folder ' + str(i+1) + ' of ' + str(len(img_dirs)) + '\n')
                sys.stderr.flush()                
            ko_fp = os.path.join(args.IMG_dir, imgID, imgID + ".ko.tab.txt")
            if os.path.exists(ko_fp):
                with open(ko_fp, 'rb') as inf:
                    database_gen = read_img_ko_table(inf, imgID)
                    for gene_id, ko in database_gen:
                        outf.write(gene_id + '\t' + ko + '\n')

if __name__ == '__main__':
    main()
