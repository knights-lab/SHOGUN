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


def intersection():
    pass


def main():
    parser = make_arg_parser()
    args = parser.parse_args()

    img2kegg = defaultdict(set)
    with open(args.mapping) as inf:
        csv_inf = csv.reader(inf, delimiter='\t')
        for row in csv_inf:
            img2kegg[row[0][:row[0].find('_')]].add(row[1])

    # img_df = pd.read_csv(args.mapping, delimiter='\t', header=None, index_col=0, names=['img_id', 'kegg_id'])
    # img_df['img_species_id'] = [i[:i.find('_')] for i in img_df.index.values]
    # for index, row in img_df.iterrows():
    #     img2kegg[row['img_species_id']].add(row['kegg_id'])

    counts = defaultdict(lambda: defaultdict(int))

    with open(args.input) as inf:
        csv_inf = csv.reader(inf)
        header = next(csv_inf)
        # 'sample_id', 'sequence_id', 'ncbi_tid', 'img_id'
        for line in csv_inf:
            img_ids = line[-1].split(',')
            sample_dict = counts[line[0]]
            if img_ids:
                if len(img_ids) > 1:
                    first = img_ids.pop()
                    kegg_intersection = img2kegg[first[:first.find('_')]]
                    for img_id in img_ids:
                        kegg_ids = img2kegg[img_id[:img_id.find('_')]]
                        kegg_intersection = kegg_ids.intersection(kegg_intersection)
                else:
                    kegg_intersection = img2kegg[img_ids[0][:img_ids[0].find('_')]]
                for key in kegg_intersection:
                    sample_dict[key] += 1

    with open(args.output, 'wb') if args.output else sys.stdout as outf:
        count_df = pd.DataFrame(counts)
        outf.write(count_df.T.to_csv())



if __name__ == '__main__':
    main()
