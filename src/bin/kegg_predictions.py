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
    parser.add_argument('-a', '--algorithm', help='If nothing is given, then intersection, else stdout.',
                        choices=['intersection', 'consensus'], default='intersection')
    return parser


def split_iter(string, sep=''):
    r = re.compile(sep)
    return (x.group(0) for x in re.finditer(r, string))


def intersection(csv_imgs_ids_inf, img_oid2kegg):
    """
    Returns the intersection of kegg_ids for each sequence for each sample in the study

    :param csv_imgs_ids_inf: the file iter from the kegg_parse_img_ids script
    :param img_oid2kegg: the img2kegg mapping dictionary
    :return:
    """

    counts = defaultdict(lambda: defaultdict(int))
    # 'sample_id', 'sequence_id', 'ncbi_tid', 'img_id'
    for line in csv_imgs_ids_inf:
        img_ids = line[-1].split(',')
        sample_dict = counts[line[0]]
        if img_ids:
            if len(img_ids) > 1:
                first = img_ids.pop()
                kegg_intersection = img_oid2kegg[first[:first.find('_')]]
                for img_id in img_ids:
                    kegg_ids = img_oid2kegg[img_id[:img_id.find('_')]]
                    kegg_intersection = kegg_ids.intersection(kegg_intersection)
            else:
                kegg_intersection = img_oid2kegg[img_ids[0][:img_ids[0].find('_')]]
            for key in kegg_intersection:
                sample_dict[key] += 1
    return counts


def get_img_oid2kegg(mapping_path):
    img_oid2kegg = defaultdict(set)
    with open(mapping_path) as inf:
        csv_inf = csv.reader(inf, delimiter='\t')
        for row in csv_inf:
            img_oid2kegg[row[0][:row[0].find('_')]].add(row[1])
    return img_oid2kegg


def get_img2kegg(mapping_path):
    with open(mapping_path) as inf:
        csv_inf = csv.reader(inf, delimiter='\t')
        img2kegg = defaultdict(str, csv_inf)
    return img2kegg


def consensus(csv_img_ids_inf, img2kegg):
    counts = defaultdict(lambda: defaultdict(int))

    # 'sample_id', 'sequence_id', 'ncbi_tid', 'img_id'
    for line in csv_img_ids_inf:
        img_ids = line[-1].split(',')
        sample_dict = counts[line[0]]
        if img_ids:
            if len(img_ids) > 1:
                first = img_ids.pop()
                kegg = img2kegg[first]
                for img_id in img_ids:
                    if kegg != img2kegg[img_id]:
                        kegg = None
                        break
            else:
                kegg = img2kegg[img_ids[0]]
            if kegg:
                sample_dict[kegg] += 1
    return counts


def main():
    parser = make_arg_parser()
    args = parser.parse_args()

    with open(args.input) as inf:
        csv_inf = csv.reader(inf)
        header = next(csv_inf)
        # 'sample_id', 'sequence_id', 'ncbi_tid', 'img_id'
        if args.algorithm in ['intersection']:
            img_oid2kegg = get_img_oid2kegg(args.mapping)
            if args.algorithm == 'intersection':
                counts = intersection(csv_inf, img_oid2kegg)
        else:
            img2kegg = get_img2kegg(args.mapping)
            if args.algorithm == 'consensus':
                counts = consensus(csv_inf, img2kegg)

    with open(args.output, 'w') if args.output else sys.stdout as outf:
        count_df = pd.DataFrame(counts)
        outf.write(count_df.T.to_csv())


if __name__ == '__main__':
    main()
