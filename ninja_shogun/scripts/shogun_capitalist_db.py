#!/usr/bin/env python
import click
import os

from ninja_utils.parsers import FASTA
from ninja_utils.utils import find_between
from ninja_utils.utils import verify_make_dir

from ninja_dojo.taxonomy import NCBITree


@click.command()
@click.option('-i', '--input', type=click.File(), default='-')
@click.option('-o', '--output', type=click.Path(), default=os.path.join(os.getcwd(), 'annotated'))
@click.option('-x', '--extract_ncbi_tid', default='ncbi_tid|,|')
def shogun_capitalist_db(input, output, extract_ncbi_tid):

    verify_make_dir(output)
    verify_make_dir(os.path.join(output, 'capitalist'))

    begin, end = extract_ncbi_tid.split(',')
    tree = NCBITree()

    inf_fasta = FASTA(input)
    for title, seq in inf_fasta.read():
        title = '>' + title
        title.replace('\t', '|')
        ncbi_tid = int(find_between(title, begin, end))
        if ncbi_tid:
            gg = tree.gg_lineage(ncbi_tid[0])
            with open(os.path.join(output, 'capitalist', gg + '.fna'), 'a+') as output_fna:
                header = 'ncbi_tid|%d|%s' % (ncbi_tid, title[1:])
                output_fna.write('>%s\n%s\n' % (header, seq))


if __name__ == '__main__':
    shogun_capitalist_db()
