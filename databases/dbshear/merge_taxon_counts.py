# combines taxon counts across multiple shogun embalmer output files. 
# if files contain the same "samples" (header names), 
# *.py taxon_counts1.tsv taxon_counts2.tsv ... outputfile.tsv
import sys
from collections import Counter, defaultdict
import pandas as pd

outf = sys.argv[-1]
taxon_fps = sys.argv[1:-1]


# counts of taxa in each sample
# Counter for each sample
counts = defaultdict(Counter) # {sample_id:{taxon1:count, taxon2:count,}}

for fp in taxon_fps:
    print(fp)
    lines = open(fp).readlines()
    header = lines[0].strip().split('\t')
    sample_ids = header[1:]
    for line in lines[1:]:
        words = line.strip().split('\t')
        taxon = words[0]
        for i,word in enumerate(words[1:]):
            counts[sample_ids[i]][taxon] += int(word)

sample_ids = sorted(counts.keys())
counts_list = [counts[sample_id] for sample_id in sample_ids]

df = pd.DataFrame(counts_list, index=sample_ids)
df.T.to_csv(outf,
            index_label='Taxon',na_rep='0',sep='\t')
