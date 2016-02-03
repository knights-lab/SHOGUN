from __future__ import print_function
import pandas as pd
import numpy as np
from collections import defaultdict

def parse_line(line):
    line = np.array(line.split('\t')[:-1])
    names = line[2].split()
    if names[1] in ('cf.', 'sp.'):
        species = ' '.join(names[1:])
    else:
        species = names[1]
    return np.append(np.delete(line, 2), [names[0], species])

def parse_taxonomy(path):
    with open(path, 'r') as f:
        df = pd.DataFrame([parse_line(line) for line in f], columns=['species_id', 'strain_id', 'domain', 'genus', 'species'])
    return df

def return_taxa_map(path):
    with open(path, 'r') as f:
        df = pd.DataFrame([parse_line(line) for line in f], columns=['species_id', 'strain_id', 'domain', 'genus', 'species'])
    m = defaultdict()
    for species_id in df['species_id']:
        m[species_id] = tuple(df[df['species_id'] == species_id][['domain', 'genus', 'species']].values[0])
    return m

def main():
    pass

if __name__ == '__main__':
    main()
