#!/usr/bin/env python
import click
import sys
import os

from ninja_utils.parsers import FASTA
from ninja_utils.utils import find_between
from ninja_utils.utils import verify_make_dir

from ninja_dojo.database import RefSeqDatabase
from ninja_dojo.taxonomy import NCBITree

@click.command()
@click.option('-i', '--input', type=click.Path(), default='-')
@click.option('-o', '--output', type=click.Path(), default=os.path.join(os.getcwd(), 'annotated'))
@click.option('-x', '--extract_refseq_id', default='ref|,|')
@click.option('--prefixes', default='*', help="Supply a comma-seperated list where the options are choices"
                                              " in ('AC', 'NC', 'NG', 'NM', 'NT', 'NW', 'NZ') e.g. NC,AC default=all")
def utree_gg_annotate(input, output, extract_refseq_id, prefixes):
    db = RefSeqDatabase()

    verify_make_dir(output)
    # check for the glob prefix
    prefixes = prefixes.split(',')

    begin, end = extract_refseq_id.split(',')
    tree = NCBITree()

    if '*' in prefixes:
        prefix_set = set([_ for _ in db.refseq_prefix_mapper.keys()])
    else:
        prefix_set = set([_ for _ in prefixes])

    output_fn = '.'.join(str(os.path.basename(input)).split('.')[:-1])
    with open(input, 'r') if input != '-' else sys.stdin as inf:
        with open(os.path.join(output, output_fn + '.annotated.fna'), 'w') as output_fna:
            with open(os.path.join(output, output_fn + '.annotated.map'), 'w') as output_map:
                inf_fasta = FASTA(inf)
                for title, seq in inf_fasta.read():
                    title = '>' + title
                    title.replace('\t', '|')
                    refseq_accession_version = find_between(title, begin, end)
                    if refseq_accession_version[:2] in prefix_set:
                        ncbi_tid = db.get_ncbi_tid_from_refseq_accession_version(refseq_accession_version)
                        if ncbi_tid:
                            gg = tree.gg_lineage(ncbi_tid[0])
                            header = 'ncbi_tid|%d|%s' % (ncbi_tid[0], title[1:])
                            output_fna.write('>%s\n%s\n' % (header, seq))
                            output_map.write('%s\t%s\n' % (header, gg))

if __name__ == '__main__':
    utree_gg_annotate()
