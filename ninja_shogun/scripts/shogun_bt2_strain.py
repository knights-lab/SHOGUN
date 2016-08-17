#!/usr/bin/env python
import click
import os
from cytoolz import valmap, map

from ninja_utils.utils import find_between,verify_make_dir

from ninja_dojo.taxonomy import NCBITree

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
@click.option('-x', '--extract_ncbi_tid', default='ncbi_tid|,|')
@click.option('-d', '--depth', type=click.INT, default=7, help='The depth of the search (7=species default, 0=No Collapse)')
@click.option('-p', '--threads', type=click.INT, default=1)
def shogun_bt2_strain(input, output, bt2_indx, extract_ncbi_tid, depth, threads):
    verify_make_dir(output)

    fna_files = [os.path.join(input, filename) for filename in os.listdir(input) if filename.endswith('.fna')]

    for fna_file in fna_files:
        sam_outf = os.path.join(output, '.'.join(str(os.path.basename(fna_file)).split('.')[:-1]) + '.sam')
        print(bowtie2(fna_file, sam_outf, bt2_indx, num_threads=threads))

    tree = NCBITree()
    begin, end = extract_ncbi_tid.split(',')

    sam_files = [os.path.join(output, filename) for filename in os.listdir(output) if filename.endswith('.sam')]
    with open(os.path.join(output, 'taxon_map.tsv'), 'w') as taxon_out:
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
            print(lca_map.items())
            # filter out null values
            map(lambda k, v: taxon_out.write('%s\t%s\n' % (k, v)), lca_map.items())


if __name__ == '__main__':
    shogun_bt2_strain()
