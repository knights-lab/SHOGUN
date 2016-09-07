#!/usr/bin/env python
import click
import os

from ninja_utils.utils import verify_make_dir
from ninja_utils.parsers import FASTA

from ninja_dojo.annotaters.refseq import refseq_annotater

from ninja_shogun.wrappers import bowtie2_build

@click.command()
@click.option('-i', '--input', type=click.Path(), default='-', help='The input FASTA file for annotating with NCBI TID (default=stdin)')
@click.option('-o', '--output', type=click.Path(), default=os.path.join(os.getcwd(), 'annotated'), help='The directory to output the formatted DB and BT2 db (default=annotated)')
@click.option('-x', '--extract_refseq_id', default='ref|,|', help='Characters that sandwich the RefSeq Accession Version in the reference FASTA (default="ref|,|")')
@click.option('--prefixes', default='*', help="Supply a comma-seperated list where the options are choices"
                                              " in ('AC', 'NC', 'NG', 'NM', 'NT', 'NW', 'NZ') e.g. NC,AC default=all")
@click.option('-d', '--depth', default=7, help="The depth to annotate the map")
@click.option('-f', '--depth-force', default=True, help="Force the depth criterion if missing annotation")
def shogun_bt2_db(input, output, extract_refseq_id, prefixes, depth, depth_force):
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
                with open(input) as inf:
                    inf_fasta = FASTA(inf)
                    annotater = refseq_annotater(inf_fasta.read(), prefixes, extract_refseq_id, depth=depth, depth_force=depth_force)
                    for lines_fna, lines_map in annotater:
                        output_fna.write(lines_fna)
                        output_map.write(lines_map)
    else:
        print("Found the output files \"%s\" and \"%s\". Skipping the annotation phase for this file." % (outf_fasta, outf_map))

    # Build the output BT2 database
    verify_make_dir(os.path.join(output, 'bt2'))
    print(bowtie2_build(outf_fasta, os.path.join(output, 'bt2', output_fn)))

if __name__ == '__main__':
    shogun_bt2_db()
