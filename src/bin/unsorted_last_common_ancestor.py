#!/usr/bin/env python
from __future__ import print_function
import argparse
from collections import Counter
import sys
from itertools import repeat
import csv


# The arg parser for this wrapper
def make_arg_parser():
    parser = argparse.ArgumentParser(description='Get least common ancestor for alignments in unsorted BAM/SAM file')
    parser.add_argument('input', type=str, help='The no header SAM file to parse.')
    parser.add_argument('-x', '--index', help='The taxonomy index.', required=True)
    parser.add_argument('-o', '--output', help='If nothing is given, then stdout, else write to file')
    parser.add_argument('-v', '--verbose', help='Print extra statistics', action='store_true', default=False)
    return parser


def get_correct_species(taxonomy):
    names = taxonomy.split()
    if len(names) > 1:
        if names[1] in ('cf.', 'sp.'):
            species = ' '.join(names[1:])
        else:
            species = names[1]
    else:
        species = 'NA'
    return names[0], species


def return_taxa_map(path):
    with open(path, 'r') as f:
        csv_file = csv.reader(f, delimiter='\t')
        # columns = ['taxon_oid', 'ncbi_taxon_id', 'domain', 'genus', 'species']
        m = {}
        r = next(csv_file)

        def parse_row(arr):
            genus, species = get_correct_species(arr[2])
            m[arr[0]] = tuple([arr[3], genus, species])

        # Checking for a header file
        if r[-1] == 'Finished':
            parse_row(r)

        for row in csv_file:
            parse_row(row)

        return m


def yield_alignments_from_sam_inf(inf):
    for i in inf:
        line = i.split('\t')
        yield line[0], line[2]


def build_lca_map(align_gen, taxa_path, verbose=False):
    if verbose:
        print("Loading taxon index.")

    taxa_map = return_taxa_map(taxa_path)

    def get_taxa_name(rname):
        # taxon_oid is in the taxa_map
        return ';'.join(taxa_map[rname.split('_')[0]])

    lca_map = {}
    for qname, rname in align_gen:
        if qname in lca_map:
            new_taxon = get_taxa_name(rname)
            current_qname = lca_map[qname]
            if current_qname and new_taxon:
                if current_qname != new_taxon:
                    lca_map[qname] = longest_common_ancestor(current_qname, new_taxon)
        else:
            lca_map[qname] = get_taxa_name(rname)
    # taxon count here
    taxon_counts = Counter(lca_map.values())
    #  normalize here
    return taxon_counts


# Given two Tuples that are taxonomic hierarchies, return their LCA
def longest_common_ancestor(taxa_1, taxa_2):
    taxa_1 = taxa_1.split(';')
    taxa_2 = taxa_2.split(';')
    lca = []
    for i in zip(taxa_1, taxa_2):
        if i[0] == i[1]:
            lca.append(i[0])
        else:
            break
    if lca:
        return ';'.join(lca)


def write_taxon_counts(taxon_counts, outf, verbose=False):
    # print results
    if verbose:
        print("Writing results...")

    columns = ['domain', 'genus', 'species', 'count']

    def map_counts(taxon):
        if taxon:
            taxa = taxon.split(';')
            return taxa + list(repeat(None, 3 - len(taxa))) + [taxon_counts[taxon]]

    wr = csv.writer(outf, quoting=csv.QUOTE_ALL)
    wr.writerow(columns)

    map(wr.writerow, [map_counts(taxon) for taxon in taxon_counts if taxon_counts[taxon] > 0])
    outf.close()


def main():
    parser = make_arg_parser()
    args = parser.parse_args()

    # parse command line

    with open(args.input, 'r') if args.input != '-' else sys.stdin as inf:
        align_gen = yield_alignments_from_sam_inf(inf)
        taxon_counts = build_lca_map(align_gen, args.index, verbose=args.verbose)
        with open(args.output, 'w') if args.output else sys.stdout as outf:
            return write_taxon_counts(taxon_counts, outf, verbose=args.verbose)

if __name__ == '__main__':
    main()
