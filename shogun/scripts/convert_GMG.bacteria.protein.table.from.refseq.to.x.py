# usage 
# *.py refseq_orgID_to_refseq_proteinID_table.txt new_mapping_refseq_protein_ID_to_other.txt output.txt
#
# file1:
# NC_000117.1     NP_219502.1     NP_219503.1     NP_219504.1     NP_219505.1
# 
# file2:
# NP_000007.1     K00249
#

import sys
from collections import defaultdict

rs2other = {} # {refseq_protein_id:other_protein_id}
rstable = defaultdict(list) # {refseq_org_id:[other_protein_id1, other_protein_id2,...]}

# read in rs2other mapping
count = 0
for line in open(sys.argv[2],'r'):
    count += 1
    if count % 100000 == 0:
        sys.stdout.write(str(count) + ' ')
        sys.stdout.flush()
    words = line.strip().split()
    rs2other[words[0]] = words[1]

sys.stdout.write('\n')

# construct full rs2other table
for line in open(sys.argv[1],'r'):
    words = line.strip().split()
    rs_org_id = words[0]
    for word in words[1:]:
        if word in rs2other:
            rstable[rs_org_id].append(rs2other[word])

# sorted orgIDs
orgIDs = sorted(rstable.keys())

# output new table
with open(sys.argv[3],'w') as f:
    for rs_org_id in orgIDs:
        f.write(rs_org_id + '\t' + '\t'.join(sorted(rstable[rs_org_id])) + '\n')

