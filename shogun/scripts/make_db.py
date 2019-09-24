#! /usr/bin/env python

# usage:
# make_db.py assembly_summary.txt outdir dbname
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

NCBI_TAXONOMY_LINK = "ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz"
T2GG_LINK = "https://github.com/knights-lab/BURST/raw/master/embalmlets/bin/t2gg"

if __name__ == "__main__":
    assf = sys.argv[1]
    outdir = sys.argv[2]
    dbname = sys.argv[3]
    
    # make dir
    print("Making output directory " + outdir)
    try:
        os.mkdir(outdir)
    except OSError:
        print ("Refusing to overwrite output dir %s" % outdir)
        raise
    
    # install t2gg and make the taxonid:taxonomy map (unless exists)
    if not os.path.exists('tid2gg.txt'):
        print("Extracting taxid 2 taxonomy map using taxonkit")
        os.system("cut -f 6 " + assf + " > taxids.txt")
        print('Creating taxonomy lineage file with taxonkit')
        os.system("taxonkit lineage -t taxids.txt > taxonkit_output.txt")
        os.system("parse_taxonkit_output.py taxonkit_output.txt tid2gg.txt")
        #to_remove = ['names.dmp','nodes.dmp','tid2gg.bin','taxdump.tar.gz']
    else:
        print("tid2gg.txt found, skipping creation")
        
    ncbi2tax = {} # ncbi taxon ID:taxonomy (output of t2gg)
    acc2tax = {} # refseq accession:taxonomy
    ftplinks = {} # refseq accession:ftp link

    # load mapping from ncbi tax ID to taxonomy
    print("Loading ncbiID:taxonomy mapping")
    with open('tid2gg.txt','r') as f:
        for line in f:
            words = line.strip().split('\t')
            ncbi2tax[words[0]] = words[1]
            
    # loop through assembly file, gather ftp link and taxonomy
    print("Loading ftp links and taxonomies for assemblies")
    with open(assf,'r') as f:
        for line in f:
            words = line.strip().split('\t')
            acc = words[0]
            taxid = words[5]
            ftp = words[19]
            acc2tax[acc] = ncbi2tax[taxid]
            ftplinks[acc] = ftp

    # write taxonomy map
    print('Writing taxonomy map')
    with open(os.path.join(outdir,dbname + '.tax'),'w') as f:
        for acc in acc2tax:
            f.write(acc + '\t' + acc2tax[acc] + '\n')
    
    # for each strain, download file, rename headers with accID
    # example header: >lcl|NC_004061.1_cds_WP_011053539.1_1
    # make it pipe-delimited, e.g. >GCF_000010525.1|WP_011053539.1
    print('Downloading and processing genomes')
    with open(os.path.join(outdir,dbname + '.fna'),'w') as f:
        for acc in ftplinks:
            sys.stdout.write(acc + ' ')
            sys.stdout.flush()
            basename = os.path.basename(ftplinks[acc])
            filename = basename + '_cds_from_genomic.fna.gz'
            os.system("wget -O " + filename + " " + ftplinks[acc] + '/' + filename)
            os.system("gunzip -f " + filename)
            filename = filename[:-3] # remove .gz
            # find-and-replace headers, write to master file
            with open(filename,'r') as g:
                header = '' # will store the current header
                seq = '' # will store the current sequence
                for line in g:
                    if line.startswith('>'):
                        # new sequence; print last one (if not empty)
                        if seq != '':
                            f.write(header + '\n' + seq + '\n')
                            seq = ''
                        header = line[:line.index(' ')] # before whitespace is header
                        comments = line.strip()[line.index(' '):] # comments after whitespace
                        header = header[:header.rfind('_')] # drop "_1" at end
                        header = header[(header.index('_')+1):] # drop second half of ncbi ID
                        header = acc + '_' + header[(header.index('_')+1):]
                    else:
                        seq += line.strip()
            f.write(header + '\n' + seq + '\n') # don't forget to write the last sequence
        sys.stdout.write('\n')
        
    #    for f in to_remove:
    #        os.remove(f)
