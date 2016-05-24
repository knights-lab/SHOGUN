import click
from ninja_shogun.taxonomy.ncbi_tree import NCBITree
nt = NCBITree()
nt.name2taxon_id
import pandas as pd
df = pd.read_csv('./Even_40M_1.gold', sep='\t', header=None)
df[0]
[nt.name2taxonid[i] for i in df[0]]
[nt.name2taxon_id[i] for i in df[0]]
[nt.name2taxon_id[i] for i in df[0]]
df[0]
df[0][0]
[nt.name2taxon_id[i.split()[-1].replace('_', ' ')] for i in df[0]]
[nt.name2taxon_id[i.split()[-1].replace('_', ' ')] for i in df[0]]
