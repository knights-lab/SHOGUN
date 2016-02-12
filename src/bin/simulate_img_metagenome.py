#!/usr/bin/env python
from __future__ import print_function
import argparse
import re
import sys
import csv
import os
import subprocess
import shutil
from scipy.stats import powerlaw
import numpy as np


# The arg parser for this wrapper
def make_arg_parser():
    parser = argparse.ArgumentParser(description='Simulate a metagenomic community from DWGSIM according to a powerlaw.')
    parser.add_argument('input', type=str, help='The species name file. One species name per line.')
    parser.add_argument('-x', '--index', help='The taxonomy index.', required=True)
    parser.add_argument('-d', '--database', help='The database. Sorted by taxon oid will '
                                                 'significantly increase speed.', required=True)
    parser.add_argument('-o', '--output', help='The output folder.')
    parser.add_argument('-v', '--verbose', help='Print extra statistics', action='store_true', default=False)
    parser.add_argument('-s', '--seed', help='The seed for the random number generator. (default=0)', default=0, type=int)
    parser.add_argument('-n', '--num', help='The number of reads to simulate', default=100, type=int)
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


def parse_names(inf):
    names = set()
    for line in inf:
        line = line.strip()
        line = line.split('_')
        names.add('%s %s' % (line[2], line[3]))
    return names


def cross_reference_names(name_set, path):
    with open(path, 'r') as f:
        csv_file = csv.reader(f, delimiter='\t')
        # columns = ['taxon_oid', 'ncbi_taxon_id', 'domain', 'genus', 'species']
        m = []
        r = next(csv_file)

        def parse_row(arr):
            genus, species = get_correct_species(arr[2])
            if genus + ' ' + species in name_set:
                name_set.remove(genus + ' ' + species)
                m.append(arr[:2] + [arr[3], genus, species])

        # Checking for a header file
        if r[-1] == 'Finished':
            parse_row(r)

        for row in csv_file:
            parse_row(row)

    return m


def read_fasta(f):
    """
    :param f: the FASTA file
    :return: tuples of (title, seq)
    """
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


def write_genomes(taxon_oids_set, database_gen, outfolder, make_genomes):
    find = re.compile(r"^[^_]*")

    taxon_oids_set_used = set()

    for title, line in database_gen:
        taxon_oid = re.search(find, title).group(0)
        taxon_old = taxon_oid

        if taxon_oid in taxon_oids_set:
            if make_genomes:
                with open(os.path.join(outfolder, taxon_oid + '.fa'), 'a') as outf:
                    # We can leverage runs to write to file faster
                    while taxon_oid == taxon_old and title:
                        if len(line) > 500:
                            outf.write('>' + title + os.linesep)
                            outf.write(line + os.linesep)
                        taxon_old = taxon_oid
                        try:
                            title, line = next(database_gen)
                            taxon_oid = re.search(find, title).group(0)
                        except StopIteration:
                            break
            taxon_oids_set_used.add(taxon_old)

    return taxon_oids_set_used


def mason(fa_infile, fq_outfile, number, prefix, seed):
    process = subprocess.Popen(['mason_simulator', '--quiet', '-n', str(number), '-ir', fa_infile, '-o',
                                fq_outfile, '--read-name-prefix', prefix, '--seed', str(seed)], bufsize=1, stdin=subprocess.PIPE, stdout=subprocess.PIPE, close_fds=True)
    process.communicate()


def dwgsim(fa_infile, fq_outfile, number, seed):
    process = subprocess.Popen(['dwgsim', '-e', '0.001', '-r', '0.001', '-q', 'f', '-c', '0', '-2', '0', '-N', str(number), '-y', '0.0', '-z', str(seed), fa_infile, fq_outfile])
    process.communicate()

def write_taxon_counts(m, counts, outf):
    columns = ['taxon_oid', 'ncbi_taxon_id',    'domain', 'genus', 'species', 'count']

    wr = csv.writer(outf, quoting=csv.QUOTE_ALL)
    wr.writerow(columns)

    map(wr.writerow, [row + [counts[i]] for i, row in enumerate(m)])
    outf.close()


def main():
    parser = make_arg_parser()
    args = parser.parse_args()

    if not os.path.exists(args.output):
        os.makedirs(args.output)

    outfolder_genomes = os.path.join(args.output, 'genomes')

    make_genomes = True
    if not os.path.exists(outfolder_genomes):
        os.makedirs(outfolder_genomes)
    else:
        make_genomes = False

    outfolder_simulations = os.path.join(args.output, 'simulations')

    if not os.path.exists(outfolder_simulations):
        os.makedirs(outfolder_simulations)
    else:
        shutil.rmtree(outfolder_simulations)
        os.makedirs(outfolder_simulations)

    outfile_glob = os.path.join(args.output, 'simulated.fa')

    # parse command line
    with open(args.input, 'rb') if args.input else sys.stdin as inf:
        name_set = parse_names(inf)

        m = cross_reference_names(name_set, args.index)
        print(m)
        taxon_oids_set = set([row[0] for row in m])

        with open(args.database, 'rb') as d:
            taxon_oids_set_used = write_genomes(taxon_oids_set, read_fasta(d), outfolder_genomes, make_genomes)
        taxon_oids_used = list(taxon_oids_set_used)

    a = .5
    np.random.seed(args.seed)
    size = len(taxon_oids_used)
    r = powerlaw.rvs(a, size=size)
    r = (r/sum(r))*args.num
    r = np.ceil(r).astype(np.int32)

    for i, filename in enumerate(taxon_oids_used):
        dwgsim(os.path.join(outfolder_genomes, filename + '.fa'), os.path.join(outfolder_simulations, filename), r[i], args.seed)

    with open(outfile_glob, 'wb') as outfile:
        for filename in taxon_oids_used:
            filename = os.path.join(outfolder_simulations, filename + '.bwa.read1.fastq')
            with open(filename, 'rb') as readfile:
                shutil.copyfileobj(readfile, outfile)

    taxon_oids_set_used = dict(zip(taxon_oids_used, range(len(taxon_oids_used))))
    m_2 = [None]*len(taxon_oids_used)
    for i, row in enumerate(m):
        if row[0] in taxon_oids_set_used:
            m_2[taxon_oids_set_used[row[0]]] = row

    with open(os.path.join(args.output, 'expected_counts.csv'), 'wb') as outf_counts:
        write_taxon_counts(m_2, r, outf_counts)

if __name__ == '__main__':
    main()
