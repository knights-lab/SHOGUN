#!/usr/bin/env python
import click
import os
import pyfaidx
from cytoolz import valmap
import csv
from collections import defaultdict

from ninja_utils.utils import find_between, reverse_collision_dict

from ninja_dojo.taxonomy import NCBITree

from ninja_utils.parsers import FASTA
from ninja_utils.utils import verify_make_dir

from ninja_shogun.aligners.bowtie import bowtie2


def yield_alignments_from_sam_inf(inf):
    with open(inf) as fh:
        for i in fh:
            line = i.split('\t')
            # this function yields qname, rname
            try:
                yield line[0], line[2]
            except IndexError:
                print('Incorrect SAM input %s' % (inf))

@click.command()
@click.option('-i', '--input', type=click.Path(), default=os.getcwd())
@click.option('-o', '--output', type=click.Path(), default=os.getcwd())
@click.option('-b', '--bt2_indx')
@click.option('-r', '--reference_fasta')
@click.option('-x', '--extract_ncbi_tid', default='ncbi_tid|,|')
@click.option('-d', '--depth', type=click.INT, default=7, help='The depth of the search (7=species default, 0=No Collapse)')
@click.option('-p', '--threads', type=click.INT, default=1)
def shogun_capitalist(input, output, bt2_indx, reference_fasta, extract_ncbi_tid, depth, threads):
    verify_make_dir(output)

    fna_files = [os.path.join(input, filename) for filename in os.listdir(input) if filename.endswith('.fna')]

    for fna_file in fna_files:
        sam_outf = os.path.join(output, '.'.join(str(os.path.basename(fna_file)).split('.')[:-1]) + '.sam')
        print(bowtie2(fna_file, sam_outf, bt2_indx, num_threads=threads))

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
    print(list(fna_faidx.keys()))

    reference_map = defaultdict(list)
    with open('.'.join(os.path.basename(reference_fasta).split('.')[:-1]) + '.map') as inf:
        tsv_in = csv.reader(inf, delimiter='\t')
        for line in tsv_in:
            reference_map[';'.join(line[1].split('; '))].append(line[0])

    # reverse the dict to feed into embalmer
    rf_faidx = pyfaidx.Fasta(reference_fasta)

    for key in lca_map_2.keys():
        for i in reference_map[key]:
            # print(rf_faidx[i])
            break

        for basename, headers in lca_map_2[key]:
            for header in headers:
                print(fna_faidx[basename][header][:])


    # Prepare for capitalist
    # verify_make_dir(os.path.join(output, 'queries'))
    # for fna_file in fna_files:
    #     with open(fna_file) as fna_fh:
    #         fna_iter = FASTA(fna_fh)
    #         for header, seq in fna_iter.read():
    #             title = header.split()[0]
    #             if title in read_map:
    #                 with open(os.path.join(output, 'queries', read_map[title] + '.fna'), 'a+') as output_fna:
    #                     output_fna.write('>%s\n%s\n' % (header, seq))


if __name__ == '__main__':
    shogun_capitalist()
