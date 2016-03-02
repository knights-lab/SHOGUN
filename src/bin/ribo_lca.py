#!/usr/bin/env python
from __future__ import print_function, division
import argparse
import sys
from collections import Counter
import csv
from itertools import repeat

from shogun.taxonomy.algorithms.last_common_ancestor import LCA
from shogun.taxonomy.ncbi_maps.img_map import IMGMap
from shogun.taxonomy.ncbi_maps.silva_map import SilvaMap
from shogun.taxonomy.ncbi_tree import NCBITree


def make_arg_parser():
    parser = argparse.ArgumentParser(description='Get least common ancestor for alignments in unsorted BAM/SAM file')
    parser.add_argument('input', type=str, help='The no header SAM file to parse.')
    parser.add_argument('-p', '--parser', help='Options are between silva and img', default='img', type=str)
    parser.add_argument('-o', '--output', help='If nothing is given, then STDOUT, else write to file')
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


def build_lca_map(align_gen, lca, rname_parse_func):
    lca_map = {}
    for qname, rname in align_gen:
        if qname in lca_map:
            new_taxon = rname_parse_func(rname)
            current_rname = lca_map[qname]
            if current_rname and new_taxon:
                if current_rname != new_taxon:
                    lca_map[qname] = lca(current_rname, new_taxon)
        else:
            lca_map[qname] = rname_parse_func(rname)
    # taxon count here
    taxon_counts = Counter(lca_map.values())
    #  normalize here
    return taxon_counts


def write_taxon_counts(taxon_counts, tree, outf,
                       ranks=['superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']):
    ranks_set = set(ranks)

    def map_counts(ncbi_taxon_id):
        taxa = tree.get_lineage(ncbi_taxon_id, ranks=ranks_set)
        return [cols[0] for cols in taxa] + list(repeat(None, len(ranks) - len(taxa))) + [taxon_counts[ncbi_taxon_id]]

    wr = csv.writer(outf, quoting=csv.QUOTE_ALL)
    wr.writerow(ranks + ['count'])

    for row in [map_counts(taxon) for taxon in taxon_counts if taxon_counts[taxon] > 0 and taxon]:
        wr.writerow(row)


def main():
    parser = make_arg_parser()
    args = parser.parse_args()

    with open(args.input, 'r') if args.input != '-' else sys.stdin as inf:
        align_gen = yield_alignments_from_sam_inf(inf)
        if args.parser == 'silva':
            silva_map = SilvaMap.load()

            def rname_parse_func(x):
                return rname_silva_parser(x, silva_map)
        else:

            img_map = IMGMap.load()

            def rname_parse_func(x):
                return rname_img_parser(x, img_map)
        ncbi_tree = NCBITree.load()
        lca = LCA(ncbi_tree)
        taxon_counts = build_lca_map(align_gen, lca, rname_parse_func)
        with open(args.output, 'w') if args.output else sys.stdout as outf:
            write_taxon_counts(taxon_counts, ncbi_tree, outf)


if __name__ == '__main__':
    main()
