#!/usr/bin/env python
from __future__ import print_function
import argparse
import re
import csv
from collections import defaultdict
import sys


def make_arg_parser():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('-x', '--index', help='The taxonomy index.', required=True)
    parser.add_argument('-d', '--database', help='The database file to parse.', required=True)
    parser.add_argument('-o', '--output', help='If nothing is given, then stdout, else write to file')
    return parser


def read_fasta(f):
    title = None
    data = None
    for line in f:
        if line[0] == ">":
            if title:
                yield (title.strip(), data)
            title = line[1:]
            data = ''
        else:
            data += line.strip()
    if not title:
        yield None
    yield (title.strip(), data)


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


def write_gene_length_counts(outf, lengths):
    columns = ['domain', 'genus', 'species', 'length']

    wr = csv.writer(outf, quoting=csv.QUOTE_ALL)
    wr.writerow(columns)

    map(wr.writerow, [list(species) + [lengths[species]] for species in lengths])
    outf.close()


def main():
    parser = make_arg_parser()
    args = parser.parse_args()

    find = re.compile(r"^[^_]*")

    taxa_map = return_taxa_map(args.index)

    lengths = defaultdict(int)

    with open(args.database, 'rb') as inf:
        database_gen = read_fasta(inf)
        for title, line in database_gen:
            taxon_oid = re.search(find, title).group(0)
            lengths[taxa_map.get(taxon_oid)] += len(line)

    with open(args.output, 'wb') if args.output else sys.stdout as outf:
        write_gene_length_counts(outf, lengths)


if __name__ == '__main__':
    main()
