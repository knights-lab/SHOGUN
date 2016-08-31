#!/usr/bin/env python
import click
import os

from ninja_utils.utils import verify_make_dir

from ninja_dojo.scripts.annotate_fasta import annotate_fasta

from ninja_shogun.wrappers import bowtie2_build

@click.command()
@click.option('-i', '--input', type=click.Path(), default='-', help='The input FASTA file for annotating with NCBI TID (default=stdin)')
@click.option('-o', '--output', type=click.Path(), default=os.path.join(os.getcwd(), 'annotated'), help='The directory to output the formatted DB and BT2 db (default=annotated)')
@click.option('-x', '--extract_refseq_id', default='ref|,|', help='Characters that sandwich the RefSeq Accession Version in the reference FASTA (default="ref|,|")')
@click.option('--prefixes', default='*', help="Supply a comma-seperated list where the options are choices"
                                              " in ('AC', 'NC', 'NG', 'NM', 'NT', 'NW', 'NZ') e.g. NC,AC default=all")
def shogun_bt2_db(input, output, extract_refseq_id, prefixes):
    # Verify the FASTA is annotated
    if input == '-':
        output_fn = 'stdin'
    else:
        output_fn = '.'.join(str(os.path.basename(input)).split('.')[:-1])

    outf_fasta = os.path.join(output, output_fn + '.annotated.fna')
    if not os.path.isfile(outf_fasta):
        annotate_fasta(input, output, extract_refseq_id, prefixes)
    else:
        print("Found the output file \"%s\". Skipping the annotation phase for this file." % outf_fasta)

    # Build the output BT2 database
    verify_make_dir(os.path.join(output, 'bt2'))
    print(bowtie2_build(outf_fasta, os.path.join(output, 'bt2', output_fn)))

if __name__ == '__main__':
    shogun_bt2_db()
