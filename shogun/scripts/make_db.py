#! /usr/bin/env python

# usage:
# make_db.py genes assembly_summary.txt outdir dbname
# or:
# make_db.py genomes assembly_summary.txt outdir dbname
#
# dependencies:
# Linux
# wget
#
# note: asssembly summary may be obtained from ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/assembly_summary.txt
#
# note: if tid2gg.txt already exists, will not create a new one

import sys
import os
from shogun.utils.refseq import make_refseq_fasta_and_taxonomy
#from shogun.utils.ontologies import create_gene_ontology

if __name__ == "__main__":
    dbtype = sys.argv[1]
    assemblypath = sys.argv[2]
    outdir = sys.argv[3]
    dbname = sys.argv[4]
    dbpath = os.path.join(outdir,dbname + '.fna') # full path
    taxpath = os.path.join(outdir,dbname + '.tax') # full path
    kotaxpath = os.path.join(outdir,dbname + '-ko.tax') # full path

    # download the raw CDS from genomes
    coding_only = dbtype == 'genes'
    make_refseq_fasta_and_taxonomy(assemblypath,dbpath,taxpath,coding_only=coding_only)

    # create the gene ontology (taxonomy format) file from sequence headers
    #create_gene_ontology(dbpath,genepath,ko2pathwaypath=None,idmappingpath=None)

    # build the db using requested tool
    # build_db(dbpath,taxpath,aligner='burst')

    # perform shearing/mapping operation
    # make_sheared_bayes_file(dbpath,shearedbayespath,aligner='burst')
    
    # create yaml file
    # make_db_yaml_file(dbpath,taxpath,genepath)
    
