# usage: python me.py \
#  alignment.burst.otu.txt db.tax sheared_bayes.txt

import os
import sys
import csv
import pandas as pd
import numpy as np
from scipy.sparse import csr_matrix

with open(sys.argv[1], 'r') as inf:
    csv_inf = csv.reader(inf, delimiter="\t")
    columns = next(csv_inf)
    columns = dict(zip(columns[1:], range(len(columns))))
    indptr = [0]
    indices = np.array([], dtype=int)
    data = np.array([], dtype=int)
    names = []
    for ix, row in enumerate(csv_inf):
        if ix % 1000 == 0:
            print(ix)
        names.append(row[0])
        np_row = np.array(row[1:], dtype=int)
        temp_indx = [np_row > 0]
        data = np.concatenate((data, np_row[temp_indx]))
        indices = np.concatenate((indices, np.where(temp_indx)[1]))
        indptr.append(indices.shape[0])

csr = csr_matrix((data, indices, indptr), dtype=int).T

with open(sys.argv[2]) as inf:
    csv_inf = csv.reader(inf, delimiter='\t')
    name2taxonomy = dict(csv_inf)

cols_tax = [name2taxonomy[name] for name in names]
rows_tax = [name2taxonomy[_.replace(".", "_", 1)] for _ in sorted(columns, key=columns.get)]

def index_lca(str1, str2):
    for i, (s1, s2) in enumerate(zip(str1.split(';'), str2.split(';'))):
        if s1 != s2:
            return i
    return 8

dat = np.zeros((len(rows_tax), 9), dtype=int)
for i, row_name in enumerate(rows_tax):
    row = csr.getrow(i)
    for j, indx in enumerate(row.indices):
        dat[i, index_lca(rows_tax[i], cols_tax[indx])] += row.data[j]

print(str(dat[:, 0].sum()))

df = pd.DataFrame(dat, index=rows_tax)
df['sum'] = dat.sum(axis=1)
df.drop(0, axis=1, inplace=True)
df.to_csv(sys.argv[3], header=False, sep='\t')

uniqueness_rate_per_level = np.zeros(8, dtype=float)
for i in range(0, 8):
    # Take the sum of those columns
    num_hits =  df.iloc[:, i].sum()
    # Total number of possible hits
    total_hits = df['sum'].sum()
    # Uniqueness Rate
    uniqueness_rate_per_level[i] = num_hits/total_hits
levels = ['kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species', 'strain']
list(zip(levels, uniqueness_rate_per_level))
print(uniqueness_rate_per_level.sum())

