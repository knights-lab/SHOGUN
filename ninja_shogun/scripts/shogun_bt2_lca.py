#!/usr/bin/env python
from __future__ import print_function, division

import click
from collections import Counter
import os
import pandas as pd

from ninja_utils.utils import find_between
from ninja_utils.utils import verify_make_dir

from ninja_dojo.taxonomy import NCBITree

from ninja_shogun.aligners.bowtie import bowtie2


def yield_alignments_from_sam_inf(inf):
    for i in inf:
        line = i.split('\t')
        # this function yields qname, rname
        yield line[0], line[2]


def lca_gg(taxonomy_a, taxonomy_b):
    taxonomy_1 = taxonomy_a.split(';')
    taxonomy_2 = taxonomy_b.split(';')
    lca = []
    for i in zip(taxonomy_1, taxonomy_2):
        if i[0] == i[1]:
            lca.append(i[0])
        else:
            break
    if lca:
        return ';'.join(lca)


def collapse(lca_map, depth):
    if depth > 0:
        for name in lca_map:
            taxonomy = lca_map[name]
            if taxonomy:
                taxonomy = taxonomy.split(';')
                if len(taxonomy) < depth:
                    lca_map[name] = None
                elif len(taxonomy) > depth:
                    lca_map[name] = ';'.join(taxonomy[:depth])
        return lca_map
    else:
        return lca_map


@click.command()
@click.option('-i', '--input', type=click.Path(), default=os.getcwd())
@click.option('-o', '--output', type=click.Path(), default=os.getcwd())
@click.option('-b', '--bt2_indx')
@click.option('-x', '--extract_ncbi_tid', default='ncbi_tid|,|')
@click.option('-d', '--depth', type=click.INT, default=7, help='The depth of the search (7=species default, 0=No Collapse)')
@click.option('-p', '--threads', type=click.INT, default=1)
def shogun_bt2_lca(input, output, bt2_indx, extract_ncbi_tid, depth, threads):
    verify_make_dir(output)

    fna_files = [os.path.join(input, filename) for filename in os.listdir(input) if filename.endswith('.fna')]

    for fna_file in fna_files:
        sam_outf = os.path.join(output, '.'.join(str(os.path.basename(fna_file)).split('.')[:-1]) + '.sam')
        print(sam_outf)
        print(fna_file)
        print(bowtie2(fna_file, sam_outf, bt2_indx, num_threads=threads))

    tree = NCBITree()
    begin, end = extract_ncbi_tid.split(',')

    counts = []
    sam_files = [os.path.join(output, filename) for filename in os.listdir(output) if filename.endswith('.sam')]

    for sam_file in sam_files:
        lca_map = {}
        for qname, rname in yield_alignments_from_sam_inf(sam_file):
            ncbi_tid = find_between(rname, begin, end)

            if qname in lca_map:
                new_taxon = tree.gg_lineage(ncbi_tid)
                current_rname = lca_map[qname]
                if current_rname and new_taxon:
                    if current_rname != new_taxon:
                        lca_map[qname] = lca_gg(current_rname, new_taxon)
            else:
                lca_map[qname] = tree.gg_lineage(ncbi_tid)

        lca_map = collapse(lca_map, depth)
        taxon_counts = Counter(filter(None, lca_map.values()))
        counts.append(taxon_counts)

    df = pd.DataFrame(counts, index=['#' + os.path.basename(sample)[:-3] for sample in sam_files])
    df.T.to_csv(os.path.join(output, 'taxon_counts.csv'))


if __name__ == '__main__':
    shogun_bt2_lca()
