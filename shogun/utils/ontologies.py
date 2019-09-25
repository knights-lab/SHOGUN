# /usr/bin/env python
#
# downloads mapping from KO to KEGG pathways, modules
#
# usage: 8
# download_and_parse_kegg_orthologies.py idmapping.dat dbpath genepath ko2pathwayfilepath  output.txt
#
# note: if idmapping.dat file does not exist, will download

import os
import sys
import gzip
from collections import defaultdict

IDMAPPING_LINK = "ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/idmapping.dat.gz"

# returns KO 2 KEGG Pathway map
# returns dict unless outfile provided.
def get_ko2pathway_map(outfile=None):
    os.system('wget http://rest.kegg.jp/get/br:ko00001 -O kegg_pathway_htext.txt')
    print('Parsing KEGG pathway htext.')
    ko2pathway = {}
    # start of file is
    # +D    KO
    # !
    # A09100 Metabolism
    # B  09101 Carbohydrate metabolism
    # C    00010 Glycolysis / Gluconeogenesis [PATH:ko00010]
    # D      K00844  HK; hexokinase [EC:2.7.1.1]
    # D      K12407  GCK; glucokinase [EC:2.7.1.2]
    # D      K00845  glk; glucokinase [EC:2.7.1.2]
    # ...
    # !
    #
    # [ KO | BRITE | KEGG2 | KEGG ]
    # Last updated: September 25, 2019
    #
    # Note: entries can sometimes have multiple ECs:
    # D      K07250  gabT; 4-aminobutyrate aminotransferase / (S)-3-amino-2-methylpropionate transaminase / 5-aminovalerate transaminase [EC:2.6.1.19 2.6.1.22 2.6.1.48]
    #
    with open('kegg_pathway_htext.txt','r') as f:
        pathway = ['','','',''] # will contain list of entries L1-L4
        for line in f:
            if line.startswith('A'):
                # e.g. A09100 Metabolism
                pathway[0] = '1. ' + ' '.join(line[1:].strip().split()[1:])
            elif line.startswith('B'):
                # e.g. B  09101 Carbohydrate metabolism
                pathway[1] = '2. ' + ' '.join(line[1:].strip().split()[1:])
            elif line.startswith('C'):
                # e.g. C    00010 Glycolysis / Gluconeogenesis [PATH:ko00010]
                pathway[2] = '3. ' + ' '.join(line[1:].strip().split()[1:])
            elif line.startswith('D'):
                # e.g. # D      K00845  glk; glucokinase [EC:2.7.1.2]
                pathway[3] = '4. ' + ' '.join(line[1:].strip().split()[1:])
                ko = pathway[3].split()[0]
                # remove quotes or apostrophes from text because
                # they can cause file parsing problems later
                pathway_string = ';'.join(pathway)
                pathway_string = pathway_string.replace("'","").replace('"','')
                ko2pathway[ko] = pathway_string

    print(str(len(ko2pathway)) + ' KOs mapped to ' + str(len(set(ko2pathway.values()))) + ' Pathways')
    if outfile is None:
        return ko2pathway 
    else:
        keys = sorted(ko2pathway.keys())
        with open(outfile,'w') as f:
            for ko in keys:
                f.write(ko + '\t' + ko2pathway[ko] + '\n')

# uses KEGG REST server to create
# mapping from module to kos in that module
def get_module2ko_mapping(dbpath,genepath,ko2pathwaypath=None,idmappingpath=None):

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


# uses KEGG REST server to create
# mapping from ko to module containing that KO
def get_ko2module_mapping(dbpath,genepath,ko2pathwaypath=None,idmappingpath=None):
    m2k = get_module2ko_mapping(dbpath,genepath)
    k2m = defaultdict(set)
    for module in m2k:
        for k in m2k[module]:
            k2m[k].add(module)
    return(k2m)

# downloads or reads in idmapping.dat.gz from UniProt
# finds a mapping from refseq2 other ontologies
def get_refseq2ko_mapping(outfile=None):
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
    count = 0
    print("Parsing idmapping.dat.gz.")
    with gzip.open('idmapping.dat.gz', 'rt') as f:
        for line in f:
            count += 1
            if count % 1000000 == 0:
                sys.stdout.write(str(count) + ' ')
                sys.stdout.write('\n')
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
    sys.stdout.write('\n')

    print(str(len(refseq2ko)) + ' refseqIDs mapped to ' + str(len(set(refseq2ko.values()))) + ' KOs')
    if outfile is None:
        return(refseq2ko)
    else:
        keys = sorted(refseq2ko.keys())
        with open(outfile,'w') as f:
            for refseq in keys:
                f.write(refseq + '\t' + refseq2ko[refseq] + '\n')

if __name__ == "__main__":

#    refseq2ko = get_refseq2ko_mapping('refseq2ko.txt')
    
#    get_refseq2kegg_pathway_ontology(dbpath='tmp/tmp.fna',genepath='')
    ko2pathway = get_ko2pathway_map(outfile='ko2pathway.txt')
