#!/usr/bin/env python
from __future__ import print_function, division

import argparse
import csv
import sys
from collections import Counter, defaultdict
import os
import pandas as pd
import csv

from shogun import SETTINGS
from shogun.taxonomy.algorithms.last_common_ancestor import LCA
from shogun.taxonomy.ncbi_maps.img_map import IMGMap
from shogun.taxonomy.ncbi_tree import NCBITree


def make_arg_parser():
    parser = argparse.ArgumentParser(description='Get least common ancestor for alignments in unsorted BAM/SAM file')
    parser.add_argument('-i', '--input', help='The folder containing the SAM files to process.', required=True, type=str)
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


def rname_img_parser(rname, img_map):
    img_id = int(rname.split('_')[0])
    ncbi_taxon_id = img_map.get(img_id)
    return ncbi_taxon_id


def build_lca_map(align_gen, lca, img_map):
    lca_map = defaultdict(lambda: [set(), None])
    for qname, rname in align_gen:
        img_id = int(rname.split('_')[0])
        if qname in lca_map:
            current_rname = lca_map[qname][1]
            new_taxon = img_map(img_id)
            if current_rname and new_taxon:
                if current_rname != new_taxon:
                    lca_map[qname][1] = lca(current_rname, new_taxon)
        else:
            lca_map[qname][1] = img_map(img_id)
        lca_map[qname][0].add(rname)
    # taxon count here
    #  normalize here
    return lca_map


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


def write_taxon_counts(taxon_counts, outf):
    wr = csv.writer(outf, quoting=csv.QUOTE_ALL)
    for taxon in [taxon for taxon in taxon_counts if taxon_counts[taxon] > 0 and taxon]:
        wr.writerow((taxon, taxon_counts[taxon]))


def main():
    parser = make_arg_parser()
    args = parser.parse_args()

    sam_files = [os.path.join(args.input, filename) for filename in os.listdir(args.input) if filename.endswith('.sam')]

    img_map = IMGMap.load()

    ncbi_tree = NCBITree.load()
    lca = LCA(ncbi_tree, args.depth)

    counts = []

    with open(args.output, 'w') if args.output else sys.stdout as outf:
        csv_outf = csv.writer(outf, quoting=csv.QUOTE_ALL)
        for file in sam_files:
            with open(file) as inf:
                lca_map = build_lca_map(yield_alignments_from_sam_inf(inf), lca, img_map)
                for key in lca_map:
                    img_ids, ncbi_tid = lca_map[key]
                    csv_outf.writerow([os.path.basename(file).split('.')[0],  key, ncbi_tid, ','.join(img_ids)])

if __name__ == '__main__':
    main()
