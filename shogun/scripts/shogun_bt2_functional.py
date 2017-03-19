#!/usr/bin/env python
"""
Copyright 2015-2017 Knights Lab, Regents of the University of Minnesota.

This software is released under the GNU Affero General Public License (AGPL) v3.0 License.
"""

import click
from collections import defaultdict
import os
import pandas as pd
from cytoolz import valmap, valfilter
import csv

from ninja_utils.utils import find_between
from ninja_utils.utils import verify_make_dir

from dojo.taxonomy import NCBITree
from dojo.taxonomy.maps import IMGMap

from shogun.wrappers import bowtie2_align
from shogun.parsers import yield_alignments_from_sam_inf


def build_img_ncbi_map(align_gen, lca, img_map):
    """
    Given a generator for SAM file, return a dictionary with QNAME as the key
    and (IMG IDs: list, LCA NCBI ID: int) as the value.
    :param align_gen:
    :param lca:
    :param img_map:
    :return:
    """
    lca_map = defaultdict(lambda: [set(), None])
    for qname, rname in align_gen:
        img_id = int(rname.split('_')[0])
        if qname in lca_map:
            current_rname = lca_map[qname][1]
            new_taxon = img_map(img_id)
            if current_rname and new_taxon:
                if current_rname != new_taxon:
                    lca_map[qname][1] = lca(current_rname, new_taxon)
        else:
            lca_map[qname][1] = img_map(img_id)
        lca_map[qname][0].add(rname)
    return lca_map




@click.command()
@click.option('-i', '--input', type=click.Path(), default=os.getcwd(),
              help='Directory containing the input FASTA files with ".fna" extensions (default=cwd)')
@click.option('-o', '--output', type=click.Path(), default=os.path.join(os.getcwd(), 'shogun_bt2_lca_out'),
              help='Output directory for the results')
@click.option('-b', '--bt2_indx', help='Path to the bowtie2 index')
@click.option('-x', '--extract_ncbi_tid', default='ncbi_tid|,|',
              help='Characters that sandwich the NCBI TID in the reference FASTA (default="ncbi_tid|,|")')
@click.option('-p', '--threads', type=click.INT, default=1, help='The number of threads to use (default=1)')
def shogun_functional(input, output, bt2_indx, extract_ncbi_tid, threads):
    verify_make_dir(output)

    basenames = [os.path.basename(filename)[:-4] for filename in os.listdir(input) if filename.endswith('.fna')]

    # Create a SAM file for each input FASTA file
    for basename in basenames:
        fna_inf = os.path.join(input, basename + '.fna')
        sam_outf = os.path.join(output, basename + '.sam')
        if os.path.isfile(sam_outf):
            print("Found the samfile \"%s\". Skipping the alignment phase for this file." % sam_outf)
        else:
            print(bowtie2_align(fna_inf, sam_outf, bt2_indx, num_threads=threads))

    img_map = IMGMap()

    for basename in basenames:
        sam_inf = os.path.join(output, basename + '.sam')
        step_outf = 'test'
        if os.path.isfile(step_outf):
            print("Found the \"%s.kegg.csv\". Skipping the LCA phase for this file." % step_outf)
        else:
            lca_map = build_img_ncbi_map(yield_alignments_from_sam_inf(sam_inf), )

    sam_files = [os.path.join(args.input, filename) for filename in os.listdir(args.input) if filename.endswith('.sam')]

    img_map = IMGMap()

    ncbi_tree = NCBITree()
    lca = LCA(ncbi_tree, args.depth)

    with open(args.output, 'w') if args.output else sys.stdout as outf:
        csv_outf = csv.writer(outf, quoting=csv.QUOTE_ALL, lineterminator='\n')
        csv_outf.writerow(['sample_id', 'sequence_id', 'ncbi_tid', 'img_id'])
        for file in sam_files:
            with open(file) as inf:
                lca_map = build_lca_map(yield_alignments_from_sam_inf(inf), lca, img_map)
                for key in lca_map:
                    img_ids, ncbi_tid = lca_map[key]
                    csv_outf.writerow([os.path.basename(file).split('.')[0], key, ncbi_tid, ','.join(img_ids)])

    if run_lca:
        tree = NCBITree()
        rank_name = list(tree.lineage_ranks.keys())[depth - 1]
        if not rank_name:
            raise ValueError('Depth must be between 0 and 7, it was %d' % depth)

        begin, end = extract_ncbi_tid.split(',')

        counts = []
        for basename in basenames:
            sam_file = os.path.join(output, basename + '.sam')
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

            if annotate_lineage:
                lca_map = valmap(lambda x: tree.green_genes_lineage(x, depth=depth), lca_map)
                taxon_counts = Counter(filter(None, lca_map.values()))
            else:
                lca_map = valfilter(lambda x: tree.get_rank_from_taxon_id(x) == rank_name, lca_map)
                taxon_counts = Counter(filter(None, lca_map.values()))
            counts.append(taxon_counts)

        df = pd.DataFrame(counts, index=basenames)
        df.T.to_csv(os.path.join(output, 'taxon_counts.csv'))


if __name__ == '__main__':
    shogun_functional()
