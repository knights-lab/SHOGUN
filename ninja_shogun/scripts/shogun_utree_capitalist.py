#!/usr/bin/env python
import click
import os
import pyfaidx
import csv
from collections import defaultdict
import tempfile
import pandas as pd
import numpy as np

from ninja_utils.utils import find_between
from ninja_utils.utils import verify_make_dir

from ninja_shogun.wrappers import utree_search, embalmer_align


@click.command()
@click.option('-i', '--input', type=click.Path(), default=os.getcwd(), help='Directory containing the input FASTA files with ".fna" extensions (default=cwd)')
@click.option('-o', '--output', type=click.Path(), default=os.getcwd(), help='Output directory for the results')
@click.option('-u', '--utree_indx', required=True, help='Path to the bowtie2 index')
@click.option('-r', '--reference_fasta', required=True, help='Path to the annotated Reference FASTA file with ".fna" extension')
@click.option('-p', '--threads', type=click.INT, default=1, help='The number of threads to use (default=1)')
def shogun_utree_capitalist(input, output, utree_indx, reference_fasta, threads):
    verify_make_dir(output)

    fna_files = [os.path.join(input, filename) for filename in os.listdir(input) if filename.endswith('.fna')]

    for fna_file in fna_files:
        tsv_outf = os.path.join(output, '.'.join(str(os.path.basename(fna_file)).split('.')[:-1]) + '.tsv')
        print(utree_search(utree_indx, fna_file, tsv_outf))

    tsv_files = [os.path.join(output, filename) for filename in os.listdir(output) if filename.endswith('.tsv')]
    lca_maps = defaultdict(lambda: defaultdict(list))
    for tsv in tsv_files:
        basename = '.'.join(os.path.basename(tsv).split('.')[:-1])
        with open(tsv) as inf:
            tsv_parser = csv.reader(inf, delimiter='\t')
            for line in tsv_parser:
                if line[1]:
                    lca_maps[';'.join(line[1].split('; '))][basename].append(line[0])

    fna_faidx = {}
    for fna_file in fna_files:
        fna_faidx[os.path.basename(fna_file)[:-4]] = pyfaidx.Fasta(fna_file)

    reference_map = defaultdict(list)
    with open('.'.join(os.path.basename(reference_fasta).split('.')[:-1]) + '.map') as inf:
        tsv_in = csv.reader(inf, delimiter='\t')
        for line in tsv_in:
            reference_map[';'.join(line[1].split('; '))].append(line[0])

    # reverse the dict to feed into embalmer
    references_faidx = pyfaidx.Fasta(reference_fasta)

    tmpdir = tempfile.mkdtemp()
    with open(os.path.join(output, 'embalmer_out.txt'), 'w') as embalmer_cat:
        for species in lca_maps.keys():

            queries_fna_filename = os.path.join(tmpdir, 'queries.fna')
            references_fna_filename = os.path.join(tmpdir, 'reference.fna')
            output_filename = os.path.join(tmpdir, 'output.txt')

            with open(queries_fna_filename, 'w') as queries_fna:
                for basename in lca_maps[species].keys():
                    for header in lca_maps[species][basename]:
                        record = fna_faidx[basename][header][:]
                        queries_fna.write('>%s\n%s\n' % (record.name, record.seq))

            with open(references_fna_filename, 'w') as references_fna:
                for i in reference_map[species]:
                        record = references_faidx[i][:]
                        references_fna.write('>%s\n%s\n' % (record.name, record.seq))

            embalmer_align(queries_fna_filename, references_fna_filename, output_filename)

            with open(output_filename) as embalmer_out:
                for line in embalmer_out:
                    embalmer_cat.write(line)

            os.remove(queries_fna_filename)
            os.remove(references_fna_filename)
            os.remove(output_filename)

    os.rmdir(tmpdir)

    sparse_ncbi_dict = defaultdict(dict)

    # build query by NCBI_TID DataFrame
    with open(os.path.join(output, 'embalmer_out.txt')) as embalmer_cat:
        embalmer_csv = csv.reader(embalmer_cat, delimiter='\t')
        for line in embalmer_csv:
            # line[0] = qname, line[1] = rname, line[2] = %match
            ncbi_tid = np.int(find_between(line[1], begin, end))
            sparse_ncbi_dict[line[0]][ncbi_tid] = np.float(line[2])

    df = pd.DataFrame.from_dict(sparse_ncbi_dict)
    df.to_csv(os.path.join(output, 'strain_alignments.csv'))

if __name__ == '__main__':
    shogun_utree_capitalist()
