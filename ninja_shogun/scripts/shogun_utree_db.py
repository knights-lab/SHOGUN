#!/usr/bin/env python
import click
import os

from ninja_utils.utils import verify_make_dir
from ninja_utils.parsers import FASTA

from ninja_dojo.annotaters.refseq import refseq_annotater

from ninja_shogun.wrappers import utree_build, utree_compress
from ninja_shogun import SETTINGS


@click.command()
@click.option('-i', '--input', type=click.Path(), default='-', help='The input FASTA file for annotating with NCBI TID (default=stdin)')
@click.option('-o', '--output', type=click.Path(), default=os.path.join(os.getcwd(), 'annotated'), help='The directory to output the formatted DB and BT2 db (default=annotated)')
@click.option('-x', '--extract_refseq_id', default='ref|,|', help='Characters that sandwich the RefSeq Accession Version in the reference FASTA (default="ref|,|")')
@click.option('-p', '--threads', type=click.INT, default=SETTINGS.N_jobs, help='The number of threads to use (default=MAX_THREADS)')
@click.option('--prefixes', default='*', help="Supply a comma-seperated list where the options are choices"
                                              " in ('AC', 'NC', 'NG', 'NM', 'NT', 'NW', 'NZ') e.g. NC,AC default=all")
def shogun_utree_db(input, output, extract_refseq_id, threads, prefixes):
    verify_make_dir(output)
    # Verify the FASTA is annotated
    if input == '-':
        output_fn = 'stdin'
    else:
        output_fn = '.'.join(str(os.path.basename(input)).split('.')[:-1])

    outf_fasta = os.path.join(output, output_fn + '.annotated.fna')
    outf_map = os.path.join(output, output_fn + '.annotated.map')
    if not os.path.isfile(outf_fasta) or not os.path.isfile(outf_map):
        with open(outf_fasta, 'w') as output_fna:
            with open(outf_map, 'w') as output_map:
                inf_fasta = FASTA(outf_fasta)
                annotater = refseq_annotater(inf_fasta.read(), prefixes, extract_refseq_id, depth=depth, depth_force=depth_force)
                for lines_fna, lines_map in annotater:
                    output_fna.write(lines_fna)
                    output_map.write(lines_map)

    # Build the output CTR
    verify_make_dir(os.path.join(output, 'utree'))
    path_uncompressed_tree = os.path.join(output, 'utree', output_fn + '.utr')
    path_compressed_tree = os.path.join(output, 'utree', output_fn+'.ctr')
    if os.path.exists(path_compressed_tree):
        print('Compressed tree database file %s exists, skipping this step.' % path_compressed_tree)
    else:
        if not os.path.exists(path_uncompressed_tree):
            print(utree_build(outf_fasta, outf_map, path_uncompressed_tree, threads=threads))
        print(utree_compress(path_uncompressed_tree, path_compressed_tree))
        os.remove(path_uncompressed_tree)

if __name__ == '__main__':
    shogun_utree_db()
