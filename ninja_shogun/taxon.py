from __future__ import print_function
from collections import defaultdict

def load_taxon_map(path, verbose=False):
    if verbose:
        # load the taxon map
        print("Loading taxonomy map...")
    taxa = defaultdict(str)
    with open(path, 'r') as f:
        for line in f:
            words = line.strip().split('\t')
            taxonomy = tuple(words[1].split(';'))
            if len(taxonomy) == 1 and words[1] != 'Archaea' and words[1] != 'Eukaryota':
                taxonomy = tuple(taxonomy[0].split(' '))
            taxa[words[0]] = taxonomy
    return taxa

if __name__ == '__main__':
    import os
    taxon_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), '..', '..', '..', 'data', 'img.taxonomy.txt')
    print(load_taxon_map(taxon_path))
