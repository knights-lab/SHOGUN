#!/usr/bin/env python
from __future__ import print_function
import pysam
import argparse
from sheq.parse_taxonomy import return_taxa_map
from collections import defaultdict
import sys
import multiprocessing
import pandas as pd
from itertools import repeat


# The arg parser for this wrapper
def make_arg_parser():
    parser = argparse.ArgumentParser(description='Get least common ancestor for alignments in sorted BAM/SAM file')
    parser.add_argument('input', help='If nothing is given then stdin, else the input BAM/SAM file',
                        default='-')
    parser.add_argument('-x', '--index', help='The taxonomy index.', required=True)
    parser.add_argument('-o', '--output', help='If nothing is given, then stdout, else write to file')
    parser.add_argument('-v', '--verbose', help='Print extra statistics', action='store_true', default=False)
    parser.add_argument('-p', '--processors', help='The number of parallel cores to use', default=None, type=int)
    return parser

def load_sam(path):
    return pysam.AlignmentFile(path, 'r')

# Under the assumption that the BAM/SAM file is sorted by query_name
def yield_all_alignments_for_queries(sam_path, taxa_path, verbose=False):
    samfile = load_sam(sam_path)

    if verbose:
        print("Loading taxon index.")

    taxa_map = return_taxa_map(taxa_path)

    def get_taxa_name(align):
        get_rname = samfile.getrname(align.rname).split('_')
        return taxa_map[get_rname[0]]

    root = samfile.next()
    current_que = root.query_name
    current_refs = set()
    current_refs.add(get_taxa_name(root))
    for i, align in enumerate(samfile):
        if verbose and i % 100000 == 0:
            print(current_que, current_refs, i)
        if align.query_name != current_que:
            yield current_que, current_refs
            current_que = align.query_name
            current_refs = set()
        current_refs.add(get_taxa_name(align))

# given two tuples that are taxonomic hierarchies,
# returns the deepest shared hierarchy
# or None if no intersection
def least_common_ancestor(taxa_set):
    for i, level in enumerate(reversed(zip(*taxa_set))):
        if len(set(level)) == 1:
            convert_i = lambda x: None if x == 0 else -x
            return ';'.join([taxon_level[0] for taxon_level in list(zip(*taxa_set))[:convert_i(i)]])
    return None

def parse_alignments(align_tuple):
    que_name, taxa_set = align_tuple
    taxon = least_common_ancestor(taxa_set)
    return que_name, taxon

def write_taxon_counts(taxon_assignments, outf, verbose=False):
    # tabulate taxon counts
    if verbose:
        print(len(taxon_assignments), "hits processed. Tabulating taxa...")
    hits = 0
    taxon_counts = defaultdict(int)
    for tuple in taxon_assignments:
        que_name, taxon = tuple
        if taxon:
            taxon_counts[taxon] += 1
            hits += 1
    if verbose:
        print(float(hits) / len(taxon_assignments) * 100, "percent of hits matched our criteria.")

    # print results
    if verbose:
        print("Writing results...")

    columns = ['domain', 'genus', 'species', 'count']

    def map_counts(taxon):
        taxa = taxon.split(';')
        return taxa + list(repeat(None, 3-len(taxa))) + [taxon_counts[taxon]]

    df = pd.DataFrame([map_counts(taxon) for taxon in taxon_counts], columns=columns)
    outf.write(df.to_csv(index=False))
    outf.close()

def main():
    parser = make_arg_parser()
    args = parser.parse_args()

    # parse command line

    cpu_count = args.processors if args.processors else multiprocessing.cpu_count()
    pool = multiprocessing.Pool(cpu_count)
    taxon_assignments = pool.map(parse_alignments,
                                 yield_all_alignments_for_queries(args.input, args.index, verbose=args.verbose),
                                 chunksize=100000)
    with open(args.output, 'w') if args.output else sys.stdout as outf:
        return write_taxon_counts(taxon_assignments, outf, verbose=args.verbose)

if __name__ == '__main__':
    main()


