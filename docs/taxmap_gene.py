# usage: taxmap gene.py

import sys
import csv
import re

RefSeq2Uniprot = dict()
Uniprot2KO = dict()

with open(sys.argv[1], 'r') as inf:
    csv_inf = csv.reader(inf, delimiter="\t")
    for i, row in enumerate(csv_inf):
        if row[1] == "RefSeq":
            RefSeq2Uniprot[row[2]] = row[0]
        elif row[1] == "KO":
            Uniprot2KO[row[0]] = row[2]

with open(sys.argv[2], 'r') as inf:
    for line in inf:
        title_search = re.search('\[protein_id=(.*?)\]', line, re.IGNORECASE)
        output_str = line[1:] + "\t"
        if title_search:
            title = title_search.group(1)

            if title in RefSeq2Uniprot:
                uniprot_id = RefSeq2Uniprot[title]
                if uniprot_id in Uniprot2KO:
                    ko_id = Uniprot2KO[uniprot_id]
                    output_str += ko_id
        print(output_str)

