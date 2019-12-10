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
MODULE_LINK = 'http://rest.kegg.jp/get/br:ko00002'
EC_LINK = 'ftp://ftp.expasy.org/databases/enzyme/enzclass.txt'
PATHWAY_LINK = 'http://rest.kegg.jp/get/br:ko00001'

# creates a functional ontology map "taxonomy" file from fasta headers
def get_refseqfastq2ontology_map(fastafp, refseq2other, outfile=None, overwrite_existing_resources=False):

    """
    Example fasta header:
    >GCF_000005825.2|WP_012957018.1|1 [locus_tag=BPOF4_RS00005] [protein=chromosomal replication initiator protein DnaA] [protein_id=WP_012957018.1] [location=816..2168] [gbkey=CDS]
    """
    outmap = {} # observed refeqID:other ontology
    with open(fastqfp,'r') as f:
        for line in f:
            if not line[0] == '>':
                continue
            if not '[protein_id=' in line:
                continue
            refseqID = line[line.index('[protein_id=') + 13:]
            refseqID = refseqID[:refseqID.index(']')]
            if refseqID in refseq2other:
                

    # read file using gzip
    # example:
    # ...
    # A9MC22  RefSeq  WP_002965908.1

    
# downloads or reads in idmapping.dat.gz from UniProt
# finds a mapping from refseq2 other ontologies
def get_ontology2ontology_map(outfile=None,
                              ontology1=['RefSeq','KO','UniRef100','UniRef90','UniRef50'][0],
                              ontology2=['RefSeq','KO','UniRef100','UniRef90','UniRef50'][2],
                              overwrite_existing_resources=False):
    if not os.path.exists('idmapping.dat.gz') or overwrite_existing_resources:
        os.system('wget --retry-connrefused --waitretry=1 --read-timeout=20 --timeout=15 -t 0 ' + IDMAPPING_LINK + ' -O idmapping.dat.gz')
    
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
            words = line.strip().split('\t')
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
            if words[1] == ontology1:
                refseqID = words[2]
            elif words[1] == ontology2:
                koID = words[2]            
        # don't forget last one
        if koID != "" and refseqID != "":
            refseq2ko[refseqID] = koID
    sys.stdout.write('\n')

    print(str(len(refseq2ko)) + ' ' + ontology1 + ' IDs mapped to ' + str(len(set(refseq2ko.values()))) + ' ' + ontology2 + ' IDs')
    if outfile is None:
        return(refseq2ko)
    else:
        keys = sorted(refseq2ko.keys())
        with open(outfile,'w') as f:
            for refseq in keys:
                f.write(refseq + '\t' + refseq2ko[refseq] + '\n')


# downloads or reads in idmapping.dat.gz from UniProt
# finds a mapping from refseq2 other ontologies
def get_refseq2ko_map(outfile=None,overwrite_existing_resources=False):
    if not os.path.exists('idmapping.dat.gz') or overwrite_existing_resources:
        os.system('wget --retry-connrefused --waitretry=1 --read-timeout=20 --timeout=15 -t 0 ' + IDMAPPING_LINK + ' -O idmapping.dat.gz')
    
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
            words = line.strip().split('\t')
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


# downloads or reads in idmapping.dat.gz from UniProt
# finds a mapping from refseq2 other ontologies
def get_refseq2pathway_map(outfile=None, overwrite_existing_resources=False):
    refseq2ko = get_refseq2ko_map(overwrite_existing_resources=overwrite_existing_resources)
    ko2pathway = get_ko2pathway_map(overwrite_existing_resources=overwrite_existing_resources)
    refseq2pathway = {}
    
    keys = sorted(refseq2ko.keys())
    missing_kos = 0
    for r in keys:
        ko = refseq2ko[r]
        if ko in ko2pathway:
            pathway = ko2pathway[ko]
            refseq2pathway[r] = pathway
        else:
            missing_kos += 1
            
    print(str(len(refseq2pathway)) + ' refseqIDs mapped to ' + str(len(set(refseq2pathway.values()))) + ' Pathways')
    if missing_kos > 0:
        print("Warning: there were " + str(missing_kos) + ' KOs present in the UniProt DB but not in the KEGG pathway mapping.')
    if outfile is None:
        return(refseq2pathway)
    else:
        keys = sorted(refseq2pathway.keys())
        with open(outfile,'w') as f:
            for refseq in keys:
                f.write(refseq + '\t' + refseq2pathway[refseq] + '\n')

# returns KO 2 KEGG Pathway map
# returns dict of lists of pathways,
# unless outfile provided, in which case
# writes one ko:pathway mapping per line
# Note: some KOs show up multiple times.
def get_ko2pathway_map(outfile=None, skip=['Human Diseases','Not Included in Pathway or Brite','Organismal Systems'],
                       overwrite_existing_resources=False):
    if os.path.exists('kegg_pathway_htext.txt') and not overwrite_existing_resources:
        print('Warning: kegg_pathway_htext.txt exists; skipping download.')
    else:
        os.system('wget --retry-connrefused --waitretry=1 --read-timeout=20 --timeout=15 -t 0 ' + PATHWAY_LINK + ' -O kegg_pathway_htext.txt')
    print('Parsing KEGG pathway htext.')
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
    # Note: ignore these L1 categories:
    # Human Diseases
    # Not Included in Pathway or Brite
    # Organismal Systems
    ko2pathway = defaultdict(set)
    allpathways = set()
    with open('kegg_pathway_htext.txt','r') as f:
        pathway = ['','','',''] # will contain list of entries L1-L4
        for line in f:
            # remove any existing semi-colons
            line = line.replace(';','')

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
                if pathway[0][3:] not in skip:
                    # e.g. # D      K00845  glk; glucokinase [EC:2.7.1.2]
                    pathway[3] = '4. ' + ' '.join(line[1:].strip().split()[1:])
                    ko = line[1:].strip().split()[0]
                    # remove quotes or apostrophes from text because
                    # they can cause file parsing problems later
                    pathway_string = ';'.join(pathway)
                    pathway_string = pathway_string.replace("'","").replace('"','')
                    allpathways.add(pathway_string)
                    pathway_string += ';' + ko
                    ko2pathway[ko].add(pathway_string)                        

    print(str(len(ko2pathway)) + ' KOs mapped to ' + str(len(allpathways)) + ' Pathways')
    if outfile is None:
        return ko2pathway 
    else:
        keys = sorted(ko2pathway.keys())
        with open(outfile,'w') as f:
            for ko in keys:
                for p in ko2pathway[ko]:
                    f.write(ko + '\t' + p + '\n')

# uses KEGG REST server to create
# mapping from module to kos in that module
def get_module2ko_map(dbpath,genepath,ko2pathwaypath=None,idmappingpath=None,overwrite_existing_resources=False):

    # download kegg list of modules
    modulelist = []
    if not os.path.exists('kegg_module_list.txt') or overwrite_existing_resources:
        os.system('wget --retry-connrefused --waitretry=1 --read-timeout=20 --timeout=15 -t 0 ' + MODULE_LINK + ' -O kegg_module_list.txt')
    with open('kegg_module_list.txt','r') as f:
        for line in f:
            if line.startswith('D'):
                modulelist.append(line.split()[1])
    print(str(len(modulelist)) + " modules found.")

    # loop through modules to get module: KO mapping for each one
    ko2module = defaultdict(set)
    for module in modulelist:
        os.system('wget --retry-connrefused --waitretry=1 --read-timeout=20 --timeout=15 -t 0 http://rest.kegg.jp/get/' + module + ' -O singlemodule.txt')
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
def get_ko2module_map(dbpath,genepath,ko2pathwaypath=None,idmappingpath=None,overwrite_existing_resources=False):
    m2k = get_module2ko_map(dbpath,genepath,overwrite_existing_resources=overwrite_existing_resources)
    k2m = defaultdict(set)
    for module in m2k:
        for k in m2k[module]:
            k2m[k].add(module)
    return(k2m)


# KO to Enzyme Commission Number(s)
# Note: one KO can map to multiple ECs
def get_ko2ec_map(outfile=None,overwrite_existing_resources=False):
    # download kegg list of modules
    ko2p = get_ko2pathway_map(overwrite_existing_resources=overwrite_existing_resources)
    ko2ec = {} # {ko : [ec1, ec2,...]

    # for each KO, pull the EC out of the end of L4
    # must be the last item, in bracket
    # e.g. gabT; 4-aminobutyrate aminotransferase / (S)-3-amino-2-methylpropionate transaminase / 5-aminovalerate transaminase [EC:2.6.1.19 2.6.1.22 2.6.1.48]
    all_ecs = set()
    for ko in ko2p:
        for p in ko2p[ko]:
            l4 = p.strip().split(';')[3]
            if '[EC:' in l4:
                ecs = l4[(l4.rindex('[EC:')+4):l4.rindex(']')]
                ecs = ecs.split()
                # check that same number of gene descriptions are there
                ko2ec[ko] = ecs
                all_ecs = all_ecs.union(ecs)
#        elif '[EC' in ko2p[ko]:
#            print(ko2p[ko].strip().split(';'))

    print(str(len(ko2ec)) + ' kos mapped to ' + str(len(all_ecs)) + ' ECs')
    if outfile is None:
        return(ko2ec)
    else:
        keys = sorted(ko2ec.keys())
        with open(outfile,'w') as f:
            for ko in keys:
                for ec in ko2ec[ko]:
                    f.write(ko + '\t' + ec + '\n')
    

# KO to Enzyme Commission Number pathway
def get_ko2ecpathway_map(outfile=None, overwrite_existing_resources=False):
    # download list of EC modules
    ko2ec = get_ko2ec_map(overwrite_existing_resources=overwrite_existing_resources)
    eclevels = {} # 1. 1. 2.-  :  With a cytochrome as acceptor.
    if os.path.exists('ec_table_raw.txt') and not overwrite_existing_resources:
        print('Warning, ec_table_raw.txt exists, skipping.')
    else: 
        os.system('wget --retry-connrefused --waitretry=1 --read-timeout=20 --timeout=15 -t 0 ' + EC_LINK + ' -O ec_table_raw.txt')
    with open('ec_table_raw.txt','r') as f:
        # 1. -. -.-  Oxidoreductases.
        # 1. 1. -.-   Acting on the CH-OH group of donors.
        # 1. 1. 1.-    With NAD(+) or NADP(+) as acceptor.
        for line in f:
            if line[0] in '0123456789': # definitions start with a number
                ec = line[0:9].replace(' ','')
#                label = line[9:].strip()
                eclevels[ec] = line.strip()
    print(str(len(eclevels)) + " EC levels found.")


    # build ko 2 ec pathway map
    ko2ecpath = defaultdict(list) # {ko:[pathway1,pathway2,...],...}
    for ko in ko2ec:
        for ec in ko2ec[ko]:
            ecL1 = '.'.join(ec.split('.')[:1])+ '.-.-.-'
            ecL2 = '.'.join(ec.split('.')[:2])+ '.-.-'
            ecL3 = '.'.join(ec.split('.')[:3])+ '.-'
            ecpath = eclevels[ecL1]+ ';' + eclevels[ecL2] + ';' + eclevels[ecL3] + ';' + ec
            ko2ecpath[ko].append(ecpath)
            
    if outfile is None:
        return ko2ecpath
    else:
        with open(outfile,'w') as f:
            keys = sorted(ko2ecpath.keys())
            for ko in keys:
                for p in ko2ecpath[ko]:
                    f.write(ko + '\t' + p + '\n')

# main function included only for easy standalone testing purposes
if __name__ == "__main__":

    refseq2ko = get_ontology2ontology_map('refseq2ko.txt',ontology1='RefSeq',ontology2='KO')

#    refseq2pathway = get_refseq2pathway_map(outfile='refseq2pathway.txt')
#    get_ko2ec_map(outfile='ko2ec.txt')
#    get_ko2ecpathway_map(outfile='ko2ec.txt')
#    get_ko2pathway_map(outfile='ko2pathway.txt')

#    get_refseq2kegg_pathway_ontology(dbpath='tmp/tmp.fna',genepath='')
