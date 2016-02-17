#!/usr/bin/env python
from __future__ import print_function
import argparse
import sys
import csv
from collections import Counter

# The arg parser for this wrapper
def make_arg_parser():
    parser = argparse.ArgumentParser(description='Simulate a metagenomic community from DWGSIM according to a powerlaw.')
    parser.add_argument('-i', '--input', type=str, help='The species name file. One species name per line.', required=True)
    parser.add_argument('-o', '--output', type=str, help='The output file.')
    return parser

def main(args):
    with sys.stdin if args.input == "-" else open(args.input, 'rb') as inf:
        header = ['classifed', 'seq_id', 'taxon_id', 'seq_len', 'lca_mapping']
        csvreader = csv.reader(inf)
        with open(args.output, 'wb') if args.output else sys.stdout as outf:
            csvwriter = csv.writer(outf, quoting=csv.QUOTE_ALL)
            taxon_ids = [line[2] for line in csvreader if line[0] == 'C']
            counter_taxon_ids = Counter(taxon_ids)
            csvwriter.writerow(row) for row in counter_taxon_ids.items


if __name__ == '__main__':
    parser = make_arg_parser()
    args = parser.parse_args()
    main(args)
