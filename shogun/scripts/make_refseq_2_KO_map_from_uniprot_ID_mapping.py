# builds a mapping of refseq ID to KO 
#
# usage:
# *.py uniprot_idmapping.dat refseq2KO.txt
# 
# headers look like this:
# Q197F7  RefSeq  YP_654575.1
# P93212  KO      K06630

import sys

gene2rs = {}
gene2ko = {}
rs2ko = {}

count = 0
for line in open(sys.argv[1],'r'):
    count += 1
    words = line.strip().split('\t')
    up_id = words[0]
    if words[1] == "RefSeq":
        rs_id = words[2]
        gene2rs[up_id] = rs_id
        if up_id in gene2ko:
            rs2ko[rs_id] = gene2ko[up_id]
    elif words[1] == "KO":
        ko_id = words[2]
        gene2ko[up_id] = ko_id
        if up_id in gene2rs:
            rs2ko[gene2rs[up_id]] = ko_id

    if count % 100000 == 0:
        sys.stdout.write(str(count) + ' ')
        sys.stdout.flush()

sys.stdout.write('\n')

# sorted rsIDs
rsIDs = sorted(rs2ko.keys())

# output rs2ko file
with open(sys.argv[2],'w') as f:
    for rsID in rsIDs:
        f.write(rsID + '\t' + rs2ko[rsID] + '\n')


