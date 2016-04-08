#!/usr/bin/env python
from __future__ import print_function, division

import argparse
import csv
import sys
from collections import Counter
import os
import pandas as pd

from shogun import SETTINGS
from shogun.taxonomy.algorithms.last_common_ancestor import lca_mp
from shogun.taxonomy.ncbi_maps.img_map import IMGMap
from shogun.taxonomy.ncbi_maps.silva_map import SilvaMap
from shogun.taxonomy.ncbi_tree import NCBITree


def make_arg_parser():
    parser = argparse.ArgumentParser(description='Get least common ancestor for alignments in unsorted BAM/SAM file')
    parser.add_argument('-i', '--input', help='The folder containing the SAM files to process.', required=True, type=str)
    parser.add_argument('-p', '--parser', help='Options are between silva and img', default='img', type=str)
    parser.add_argument('-o', '--output', help='If nothing is given, then STDOUT, else write to file')
    parser.add_argument('-t', '--threads', help='The number of threads to use.', default=SETTINGS.N_jobs, type=int)
    parser.add_argument('-d', '--depth', help='The depth of the search (7=species default)', default=7, type=int)
    parser.add_argument('-v', '--verbose', help='Print extra statistics', action='store_true', default=False)
    return parser


def yield_alignments_from_sam_inf(inf):
    for i in inf:
        line = i.split('\t')
        # this function yields qname, rname
        yield line[0], line[2]


def rname_silva_parser(rname, silva_map):
    silva_acc = rname.split('.', 1)[0]
    ncbi_taxon_id = silva_map.get(silva_acc)
    return ncbi_taxon_id


def rname_img_parser(rname, img_map):
    img_id = int(rname.split('_')[0])
    ncbi_taxon_id = img_map.get(img_id)
    return ncbi_taxon_id


def build_lca_map(align_gen, lca, rname_parse_func, tree, depth):
    lca_map = {}
    for qname, rname in align_gen:
        if qname in lca_map:
            new_taxon = tree.mp_lineage(rname_parse_func(rname))
            current_rname = lca_map[qname]
            if current_rname and new_taxon:
                if current_rname != new_taxon:
                    lca_map[qname] = lca(current_rname, new_taxon)
        else:
            lca_map[qname] = tree.mp_lineage(rname_parse_func(rname))
    # taxon count here
    lca_map = collapse(lca_map, depth)
    taxon_counts = Counter(filter(None, lca_map.values()))
    print(len(list(taxon_counts.keys())))
    #  normalize here
    return taxon_counts


def collapse(lca_map, depth):
    for name in lca_map:
        taxonomy = lca_map[name]
        if taxonomy:
            taxonomy = taxonomy.split(';')
            if len(taxonomy) < depth:
                lca_map[name] = None
            elif len(taxonomy) > depth:
                lca_map[name] = ';'.join(taxonomy[:depth])
    return lca_map


def main():
    parser = make_arg_parser()
    args = parser.parse_args()

    sam_files = [os.path.join(args.input, filename) for filename in os.listdir(args.input) if filename.endswith('.sam')]
    if args.parser == 'silva':
        silva_map = SilvaMap.load()

        def rname_parse_func(x):
            return rname_silva_parser(x, silva_map)

    else:
        img_map = IMGMap.load()

        def rname_parse_func(x):
            return rname_img_parser(x, img_map)

    ncbi_tree = NCBITree.load()

    counts = []
    for file in sam_files:
        with open(file) as inf:
            counts.append(build_lca_map(yield_alignments_from_sam_inf(inf), lca_mp, rname_parse_func, ncbi_tree,
                                        args.depth))

    df = pd.DataFrame(counts, index=['#' + os.path.basename(sample).split('.')[0] for sample in sam_files])
    with open(args.output, 'w') if args.output else sys.stdout as outf:
        df.T.to_csv(outf)

if __name__ == '__main__':
    main()
