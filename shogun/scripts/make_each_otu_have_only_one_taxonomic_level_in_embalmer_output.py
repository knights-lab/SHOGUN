# usage
# *.py embalmer_hits.tsv embalmer_hits_trunc.tsv

import sys

infp = sys.argv[1]
outfp = sys.argv[2]

# read through once and build map refID:LCA
ref2lca = {}

for line in open(infp):
    words = line.strip().split('\t')
    refID = words[1]
    taxonomy = words[12]
    if refID not in ref2lca:
        ref2lca[refID] = taxonomy
    elif len(taxonomy) < len(ref2lca[refID]):
        ref2lca[refID] = taxonomy

# read through again and write output file
outf = open(outfp,'w')

for line in open(infp):
    words = line.strip().split('\t')
    refID = words[1]
    words[12] = ref2lca[refID]
    outf.write('\t'.join(words) + '\n')
outf.close()
