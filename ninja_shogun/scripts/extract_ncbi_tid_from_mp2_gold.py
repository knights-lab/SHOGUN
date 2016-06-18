#!/usr/bin/env python
import click
import pandas as pd
from glob import glob
import os
from collections import defaultdict

from ninja_shogun.taxonomy.ncbi.ncbi_tree import NCBITree


@click.command()
@click.argument('path', type=click.Path(exists=True))
@click.option('-v', '--verbose', is_flag=True)
def main(path, verbose):
    nt = NCBITree()
    for file in glob(os.path.join(path, '*.gold')):
        df = pd.read_csv(file, header=None, sep='\t')
        mp2_to_taxon_id = defaultdict(list)
        i = 0
        for mp2 in df[0]:
            for clade in iter(mp2.split()[::-1]):
                tid = nt.name2taxon_id[clade.replace('_', ' ')]
                if not 0 == tid:
                    mp2_to_taxon_id['metaphlan2_name'].append(mp2)
                    mp2_to_taxon_id['ncbi_taxon_id'].append(tid)
                    mp2_to_taxon_id['lineage'].append(nt.mp_lineage(tid))
                    break
                elif verbose:
                    print('%s not found' % clade)
                    i += 1

        df_out = pd.DataFrame(mp2_to_taxon_id, index=None)
        df_out.to_csv(os.path.join(path, file[:file.find('.')] + '.ncbi_map.csv'))
        if verbose:
            print('%d misses out of %d' % (i, len(mp2_to_taxon_id['metaphlan2_name'])))

if __name__ == '__main__':
    main()
