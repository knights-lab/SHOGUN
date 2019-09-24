# /usr/bin/env python
#
# downloads mapping from KO to KEGG pathways, modules
#
# usage:
# download_and_parse_kegg_orthologies.py idmapping.dat dbpath genepath ko2pathwayfilepath  output.txt
#
# note: if idmapping.dat file does not exist, will download

import os
import sys
import gzip
from collections import defaultdict

IDMAPPING_LINK = "ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/idmapping.dat.gz"

def create_gene_ontology(dbpath,genepath,ko2pathwaypath=None,idmappingpath=None):

    # download kegg pathway raw file
    os.system('wget http://rest.kegg.jp/get/br:ko00001 -O kegg_orthology_raw.txt')

    # download kegg list of modules
    modulelist = []
    os.system('wget http://rest.kegg.jp/get/br:ko00002 -O kegg_module_list.txt')
    with open('kegg_module_list.txt','r') as f:
        for line in f:
            if line.startswith('D'):
                modulelist.append(line.split()[1])
    print(str(len(modulelist)) + " modules found.")

    # loop through modules to get module: KO mapping for each one
    ko2module = defaultdict(set)
    for module in modulelist:
        os.system('wget http://rest.kegg.jp/get/' + module + ' -O singlemodule.txt')
        with open('singlemodule.txt','r') as f:
            lines = f.readlines()
            i = 0
            while not lines[i].startswith('ORTHOLOGY'):
                i += 1
            kolist = set()
            while not lines[i].startswith('CLASS'):
                kolist_i = lines[i].split()[0].split(',')
                kolist = kolist.union(set(kolist_i))
                i += 1
            for ko in kolist:
                ko2module[ko].add(module)
            
    print(str(len(ko2module)) + " KOs assigned to modules.")


# downloads or reads in idmapping.dat.gz from UniProt
# finds a mapping from refseq2 other ontologies
def get_refseq2ko_mapping():
    if not os.path.exists('idmapping.dat.gz'):
        os.system('wget ' + IDMAPPING_LINK + ' -O idmapping.dat.gz')
    
    # read file using gzip
    # example:
    # ...
    # A9MC22  RefSeq  WP_002965908.1
    # A9MC22  KEGG    bcs:BCAN_B0743
    # ...
    count = 0
    refseq2ko = {}
    uniprotID = ""
    refseqID = ""
    koID = ""
    with gzip.open('idmapping.dat.gz', 'rt') as f:
        for line in f:
            words = line.split('\t')
            if uniprotID != words[0]:
                if uniprotID != "":
                    # if not the first one, add to DB
                    # (if ko and refseq ID present)
                    if koID != "" and refseqID != "":
                        refseq2ko[refseqID] = koID
                    # reset
                    uniprotID = words[0]
                    refseqID = ""
                    koID = ""
                else:
                    uniprotID = words[0]
            if words[1] == "RefSeq":
                refseqID = words[2]
            elif words[1] == "KO":
                koID = words[2]            
        # don't forget last one
        if koID != "" and refseqID != "":
            refseq2ko[refseqID] = koID
    print(str(len(refseq2ko)) + ' refseqIDs mapped to ' + str(len(set(refseq2ko.values()))) + ' KOs')
    return(refseq2ko)

if __name__ == "__main__":

    refseq2ko = get_refseq2ko_mapping()
    create_gene_ontology(dbpath='',genepath='')
            
