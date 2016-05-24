#!/usr/bin/env python
from __future__ import print_function, division

import argparse
import sys
from collections import defaultdict
import os
import csv

from ninja_shogun import SETTINGS
from ninja_shogun.taxonomy.algorithms.last_common_ancestor import LCA
from ninja_shogun.taxonomy.ncbi.maps.img_map import IMGMap
from ninja_shogun.taxonomy.ncbi.ncbi_tree import NCBITree


def make_arg_parser():
    parser = argparse.ArgumentParser(description='Get least common ancestor for alignments in unsorted BAM/SAM file')
    parser.add_argument('-i', '--input', help='The folder containing the SAM files to process.', required=True, type=str)
    parser.add_argument('-o', '--output', help='If nothing is given, then STDOUT, else write to file')
    parser.add_argument('-t', '--threads', help='The number of threads to use.', default=SETTINGS.N_jobs, type=int)
    parser.add_argument('-d', '--depth', help='The depth of the search (7=species default)', default=7, type=int)
    parser.add_argument('-v', '--verbose', help='Print extra statistics', action='store_true', default=False)
    return parser


def yield_alignments_from_sam_inf(inf):
    csv_inf = csv.reader(inf, delimiter='\t')
    for line in csv_inf:
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
    return lca_map


def main():
    parser = make_arg_parser()
    args = parser.parse_args()

    sam_files = [os.path.join(args.input, filename) for filename in os.listdir(args.input) if filename.endswith('.sam')]

    img_map = IMGMap()

    ncbi_tree = NCBITree()
    lca = LCA(ncbi_tree, args.depth)

    with open(args.output, 'w') if args.output else sys.stdout as outf:
        csv_outf = csv.writer(outf, quoting=csv.QUOTE_ALL, lineterminator='\n')
        csv_outf.writerow(['sample_id', 'sequence_id', 'ncbi_tid', 'img_id'])
        for file in sam_files:
            with open(file) as inf:
                lca_map = build_lca_map(yield_alignments_from_sam_inf(inf), lca, img_map)
                for key in lca_map:
                    img_ids, ncbi_tid = lca_map[key]
                    csv_outf.writerow([os.path.basename(file).split('.')[0],  key, ncbi_tid, ','.join(img_ids)])

if __name__ == '__main__':
    main()
