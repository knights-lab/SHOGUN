#!/usr/bin/env python
import click
import sys
import os

from ninja_utils.parsers import FASTA
from ninja_utils.utils import find_between
from ninja_utils.utils import verify_make_dir

from ninja_dojo.database import RefSeqDatabase
from ninja_dojo.taxonomy import NCBITree

from ninja_shogun.wrappers import utree_build, utree_compress
from ninja_shogun import SETTINGS


@click.command()
@click.option('-i', '--input', type=click.Path(), default='-', help='The input FASTA file for annotating with NCBI TID (default=stdin)')
@click.option('-o', '--output', type=click.Path(), default=os.path.join(os.getcwd(), 'annotated'), help='The directory to output the formatted DB and BT2 db (default=annotated)')
@click.option('-x', '--extract_refseq_id', default='ref|,|', help='Characters that sandwich the RefSeq Accession Version in the reference FASTA (default="ref|,|")')
@click.option('-p', '--threads', type=click.INT, default=SETTINGS.N_jobs, help='The number of threads to use (default=MAX_THREADS)')
@click.option('--prefixes', default='*', help="Supply a comma-seperated list where the options are choices"
                                              " in ('AC', 'NC', 'NG', 'NM', 'NT', 'NW', 'NZ') e.g. NC,AC default=all")
def shogun_bt2_db(input, output, extract_refseq_id, threads, prefixes):
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

    if input == '-':
        output_fn = 'stdin'
    else:
        output_fn = '.'.join(str(os.path.basename(input)).split('.')[:-1])

    with open(input, 'r') if input != '-' else sys.stdin as inf:
        path_annotated_fasta = os.path.join(output, output_fn + '.annotated.fna')
        path_annotated_map = os.path.join(output, output_fn + '.annotated.map')
        with open(path_annotated_fasta, 'w') as output_fna:
            with open(path_annotated_map, 'w') as output_map:
                inf_fasta = FASTA(inf)
                for title, seq in inf_fasta.read():
                    title = '>' + title
                    title.replace('\t', '|')
                    refseq_accession_version = find_between(title, begin, end)
                    if refseq_accession_version[:2] in prefix_set:
                        ncbi_tid = db.get_ncbi_tid_from_refseq_accession_version(refseq_accession_version)
                        if ncbi_tid:
                            gg = tree.green_genes_lineage(ncbi_tid[0])
                            gg = '; '.join(gg.split(';'))
                            header = 'ncbi_tid|%d|%s' % (ncbi_tid[0], title[1:])
                            output_fna.write('>%s\n%s\n' % (header, seq))
                            output_map.write('%s\t%s\n' % (header.split()[0], gg))

    verify_make_dir(os.path.join(output, 'utree'))
    path_uncompressed_tree = os.path.join(output, 'utree', output_fn)
    print(utree_build(path_annotated_fasta, path_annotated_map, path_uncompressed_tree, threads=threads))
    print(utree_compress(path_uncompressed_tree, os.path.join(output, 'utree', output_fn+'.ctr')))
    os.remove(path_uncompressed_tree)

if __name__ == '__main__':
    shogun_bt2_db()
