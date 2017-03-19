# builds a list of protein types in refseq organisms
# each output line is refseq_organism_ID 	gene_function1	gene_function2	...
#
# also builds a list of gene IDs to protein descriptions
# each output line is gene_ID	gene_description
#
# usage:
# *.py GMG.fasta organism2gene_table.txt gene2description_table.txt
# 
# headers look like this:
# usage
# *.py GMG.fasta
# >lcl|NC_013791.2_cds_WP_012957018.1_1 [locus_tag=BPOF4_RS00005] [protein=chromosomal replication initiation protein DnaA] [protein_id=WP_012957018.1] [location=816..2168]

import sys
from collections import defaultdict

org2prot_id = defaultdict(list) # {orgID:[prot_id1, prot_id2,...]}
prot_id2desc = {} # {prot_id:prot_description}

count = 0
for line in open(sys.argv[1],'r'):
    if line.startswith('>') and line.find('_cds_') > 0:
        count += 1
        if count % 100000 == 0:
            sys.stdout.write(str(count) + ' ')
            sys.stdout.flush()
        org = line.strip().split('|')[1]
        org = org[:org.find('_cds_')]

        ix = line.find('protein_id=')
        if ix > 0:
            prot_id = line[ix + 11:]
            prot_id = prot_id.split(']')[0]
            org2prot_id[org].append(prot_id)

            ix = line.find('protein=')
            if ix > 0:
                prot_desc = line[ix + 8:]
                prot_desc = prot_desc.split(']')[0]
                prot_id2desc[prot_id] = prot_desc
            else:
                prot_id2desc[prot_id] = "Unknown"

sys.stdout.write('\n')

# sorted orgIDs
orgIDs = sorted(org2prot_id.keys())
geneIDs = sorted(prot_id2desc.keys())

# output org2proteinID file
with open(sys.argv[2],'w') as f:
    for orgID in orgIDs:
        f.write(orgID + '\t' + '\t'.join(sorted(org2prot_id[orgID])) + '\n')

# output gene2dscription file
with open(sys.argv[3],'w') as f:
    for geneID in geneIDs:
        f.write(geneID + '\t' + prot_id2desc[geneID] + '\n')


