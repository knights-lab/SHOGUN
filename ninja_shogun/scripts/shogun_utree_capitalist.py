#!/usr/bin/env python
import click
import os
import pyfaidx
import csv
from collections import defaultdict
import tempfile
import pandas as pd
import numpy as np

from ninja_utils.utils import find_between
from ninja_utils.utils import verify_make_dir

from ninja_shogun.wrappers import utree_search, embalmer_align


@click.command()
@click.option('-i', '--input', type=click.Path(), default=os.getcwd(), help='Directory containing the input FASTA files with ".fna" extensions (default=cwd)')
@click.option('-o', '--output', type=click.Path(), default=os.path.join(os.getcwd(), 'shogun_utree_capitalist_out'), help='Output directory for the results')
@click.option('-u', '--utree_indx', required=True, help='Path to the bowtie2 index')
@click.option('-r', '--reference_fasta', required=True, help='Path to the annotated Reference FASTA file with ".fna" extension')
@click.option('-m', '--reference_map', required=True, help='Path to the annotated Reference FASTA file with ".fna" extension')
@click.option('-x', '--extract_ncbi_tid', default='ncbi_tid|,|', help='Characters that sandwich the NCBI TID in the reference FASTA (default="ncbi_tid|,|")')
@click.option('-p', '--threads', type=click.INT, default=1, help='The number of threads to use (default=1)')
def shogun_utree_capitalist(input, output, utree_indx, reference_fasta, reference_map, extract_ncbi_tid, threads):
    verify_make_dir(output)

    basenames = [os.path.basename(filename)[:-4] for filename in os.listdir(input) if filename.endswith('.fna')]

    for basename in basenames:
        fna_file = os.path.join(input, basename + '.fna')
        tsv_outf = os.path.join(output, basename + '.utree.tsv')
        if not os.path.isfile(tsv_outf):
            print(utree_search(utree_indx, fna_file, tsv_outf))
        else:
            print("Found the output file \"%s\". Skipping the alignment phase for this file." % tsv_outf)

    embalmer_outf = os.path.join(output, 'embalmer_out.txt')
    # Indexing for emblalmer
    if not os.path.isfile(embalmer_outf):
        lca_maps = defaultdict(lambda: defaultdict(list))
        for basename in basenames:
            utree_tsv = os.path.join(output, basename + '.utree.tsv')
            with open(utree_tsv) as inf:
                tsv_parser = csv.reader(inf, delimiter='\t')
                for line in tsv_parser:
                    if line[1]:
                        lca_maps[';'.join(line[1].split('; '))][basename].append(line[0])


        fna_faidx = {}
        for basename in basenames:
            fna_faidx[basename] = pyfaidx.Fasta(os.path.join(input, basename + '.fna'))

        dict_reference_map = defaultdict(list)

        with open(reference_map) as inf:
            tsv_in = csv.reader(inf, delimiter='\t')
            for line in tsv_in:
                dict_reference_map[';'.join(line[1].split('; '))].append(line[0])

        # reverse the dict to feed into embalmer
        references_faidx = pyfaidx.Fasta(reference_fasta)

        tmpdir = tempfile.mkdtemp()
        print(tmpdir)
        with open(embalmer_outf, 'w') as embalmer_cat:
            for species in lca_maps.keys():

                queries_fna_filename = os.path.join(tmpdir, 'queries.fna')
                references_fna_filename = os.path.join(tmpdir, 'reference.fna')
                output_filename = os.path.join(tmpdir, 'output.txt')

                with open(queries_fna_filename, 'w') as queries_fna:
                    for basename in lca_maps[species].keys():
                        for header in lca_maps[species][basename]:
                            record = fna_faidx[basename][header][:]
                            queries_fna.write('>filename|%s|%s\n%s\n' % (basename, record.name, record.seq))

                with open(references_fna_filename, 'w') as references_fna:
                    for i in dict_reference_map[species]:
                            record = references_faidx[i][:]
                            references_fna.write('>%s\n%s\n' % (record.name, record.seq))

                embalmer_align(queries_fna_filename, references_fna_filename, output_filename)

                with open(output_filename) as embalmer_out:
                    for line in embalmer_out:
                        embalmer_cat.write(line)

                os.remove(queries_fna_filename)
                os.remove(references_fna_filename)
                os.remove(output_filename)

        os.rmdir(tmpdir)
    else:
        print("Found the output file \"%s\". Skipping the strain alignment phase for this file." % embalmer_outf)


    # Convert the results from embalmer into CSV
    sparse_ncbi_dict = defaultdict(dict)

    begin, end = extract_ncbi_tid.split(',')
    # build query by NCBI_TID DataFrame
    with open(embalmer_outf) as embalmer_cat:
        embalmer_csv = csv.reader(embalmer_cat, delimiter='\t')
        for line in embalmer_csv:
            # line[0] = qname, line[1] = rname, line[2] = %match
            ncbi_tid = np.int(find_between(line[1], begin, end))
            sparse_ncbi_dict[line[0]][ncbi_tid] = np.float(line[2])

    df = pd.DataFrame.from_dict(sparse_ncbi_dict)
    df.to_csv(os.path.join(output, 'strain_alignments.csv'))

if __name__ == '__main__':
    shogun_utree_capitalist()
