#!/usr/bin/env python
import click
import os
import pyfaidx
from cytoolz import valmap
import csv
from collections import defaultdict
import tempfile
import pandas as pd
import numpy as np

from ninja_utils.utils import find_between, reverse_collision_dict
from ninja_utils.utils import verify_make_dir

from ninja_dojo.taxonomy import NCBITree

from ninja_shogun.wrappers import bowtie2_align, embalmer_align
from ninja_shogun.parsers import yield_alignments_from_sam_inf


@click.command()
@click.option('-i', '--input', type=click.Path(), default=os.getcwd(), help='Directory containing the input FASTA files with ".fna" extensions (default=cwd)')
@click.option('-o', '--output', type=click.Path(), default=os.getcwd(), help='Output directory for the results')
@click.option('-b', '--bt2_indx', required=True, help='Path to the bowtie2 index')
@click.option('-r', '--reference_fasta', required=True, help='Path to the annotated Reference FASTA file with ".fna" extension')
@click.option('-m', '--dict_reference_map', required=True, help='Path to the annotated Reference FASTA file with ".fna" extension')
@click.option('-x', '--extract_ncbi_tid', default='ncbi_tid|,|', help='Characters that sandwich the NCBI TID in the reference FASTA (default="ncbi_tid|,|")')
@click.option('-d', '--depth', type=click.INT, default=7, help='The depth of the search (7=species default, 0=No Collapse)')
@click.option('-p', '--threads', type=click.INT, default=1, help='The number of threads to use (default=1)')
def shogun_bt2_capitalist(input, output, bt2_indx, reference_fasta, reference_map, extract_ncbi_tid, depth, threads):
    verify_make_dir(output)

    fna_files = [os.path.join(input, filename) for filename in os.listdir(input) if filename.endswith('.fna')]

    for fna_file in fna_files:
        sam_outf = os.path.join(output, '.'.join(str(os.path.basename(fna_file)).split('.')[:-1]) + '.sam')
        print(bowtie2_align(fna_file, sam_outf, bt2_indx, num_threads=threads))

    tree = NCBITree()
    begin, end = extract_ncbi_tid.split(',')

    sam_files = [os.path.join(output, filename) for filename in os.listdir(output) if filename.endswith('.sam')]
    lca_maps = {}
    for sam_file in sam_files:
        lca_map = {}
        for qname, rname in yield_alignments_from_sam_inf(sam_file):
            ncbi_tid = int(find_between(rname, begin, end))
            if qname in lca_map:
                current_ncbi_tid = lca_map[qname]
                if current_ncbi_tid:
                    if current_ncbi_tid != ncbi_tid:
                        lca_map[qname] = tree.lowest_common_ancestor(ncbi_tid, current_ncbi_tid)
            else:
                lca_map[qname] = ncbi_tid

        lca_map = valmap(lambda x: tree.green_genes_lineage(x, depth=depth), lca_map)
        # filter out null values
        lca_maps['.'.join(os.path.basename(sam_file).split('.')[:-1])] = reverse_collision_dict(lca_map)

    for basename in lca_maps.keys():
        lca_maps[basename] = valmap(lambda val: (basename, val), lca_maps[basename])

    lca_map_2 = defaultdict(list)
    for basename in lca_maps.keys():
        for key, val in lca_maps[basename].items():
            if key:
                lca_map_2[key].append(val)

    fna_faidx = {}
    for fna_file in fna_files:
        fna_faidx[os.path.basename(fna_file)[:-4]] = pyfaidx.Fasta(fna_file)

    dict_reference_map = defaultdict(list)
    with open(reference_map) as inf:
        tsv_in = csv.reader(inf, delimiter='\t')
        for line in tsv_in:
            dict_reference_map[';'.join(line[1].split('; '))].append(line[0])

    # reverse the dict to feed into embalmer
    references_faidx = pyfaidx.Fasta(reference_fasta)

    tmpdir = tempfile.mkdtemp()
    with open(os.path.join(output, 'embalmer_out.txt'), 'w') as embalmer_cat:
        for key in lca_map_2.keys():

            queries_fna_filename = os.path.join(tmpdir, 'queries.fna')
            references_fna_filename = os.path.join(tmpdir, 'reference.fna')
            output_filename = os.path.join(tmpdir, 'output.txt')

            with open(queries_fna_filename, 'w') as queries_fna:
                for basename, headers in lca_map_2[key]:
                    for header in headers:
                        record = fna_faidx[basename][header][:]
                        queries_fna.write('>filename|%s|%s\n%s\n' % (basename, record.name, record.seq))

            with open(references_fna_filename, 'w') as references_fna:
                for i in dict_reference_map[key]:
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
    shogun_bt2_capitalist()
