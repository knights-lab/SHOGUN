#!/usr/bin/env python
from __future__ import print_function
import argparse
import re
import csv
from collections import defaultdict
import sys
import pandas as pd


def make_arg_parser():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('-i', '--input', help='The input file.', required=True)
    parser.add_argument('-m', '--mapping', help='The img to KEGG mapping file.', required=True)
    parser.add_argument('-o', '--output', help='If nothing is given, then stdout, else write to file')
    return parser


def split_iter(string, sep=''):
    r = re.compile(sep)
    return (x.group(0) for x in re.finditer(r, string))


def main():
    parser = make_arg_parser()
    args = parser.parse_args()

    # img2kegg = {}
    # with open(args.mapping) as inf:
    #     csv_inf = csv.reader(inf, delimiter='\t')
    #     for line in csv_inf:
    #         img2kegg[line[0]] = line[1]

    img2kegg = pd.read_csv(args.mapping, delimiter='\t', header=None, index_col=0, names=['img_id', 'kegg_id'])['kegg_id']
    counts = {}

    with open(args.input) as inf:
        csv_inf = csv.reader(inf)
        header = next(csv_inf)
        # 'sample_id', 'sequence_id', 'ncbi_tid', 'img_id'
        for line in csv_inf:
            img_ids = split_iter(line[3], ',')
            for img_id in img_ids:
                line[:1]


if __name__ == '__main__':
    main()
