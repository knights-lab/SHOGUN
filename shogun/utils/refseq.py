#! /usr/bin/env python

# Utilities for downloading/parsing Refseq genes and genomes
import sys
import os

NCBI_TAXONOMY_LINK = "ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz"
TAXONKIT_LINK = 'https://github.com/shenwei356/taxonkit/releases/download/v0.5.0/taxonkit_linux_amd64.tar.gz'

def get_accession2taxonomy(assemblypath,save_taxonkit_output=True,outfile=None):
    """
    Makes a mapping of RefSeq accession IDs to taxonomy strings.
    
    Note: The output tax file is not appropriate for a
    gene-split database because that would require mapping
    /loci/ (e.g. GCF_000007365.1_cds_WP_011053539.1) to taxonomy
    rather than /genomes/ (e.g. GCF_000007365.1) to taxonomy.
    For gene-split database, run get_locus2taxonomy().

    Parameters
    ----------
    assemblypath : str
        The file location of the assembly summary from refseq
        containing the accessions of interest in column one and
        the NCBI taxon ID in column 6.
        
    outfile : str, optional
        Output path for tab-delimited tax file.

    Returns
    -------
    dict
        A dict mapping accession:taxonomy (or nothing if
        outfile is passed).
    """

    if os.path.exists('taxonkit_output.txt'):
        print('Warning: Taxon kit output file taxonkit_output.txt exist; skipping creation.')
    else:
        print("Downloading and extracting taxonkit")
        os.system('wget ' + TAXONKIT_LINK)
        os.system('tar xvzf taxonkit_linux_amd64.tar.gz')
        os.system('chmod a+x taxonkit')
        
        print("Downloading and extracting taxdump.tar.gz")
        os.system('wget ' + NCBI_TAXONOMY_LINK)
        os.system('tar xvzf taxdump.tar.gz nodes.dmp names.dmp delnodes.dmp merged.dmp')
        
        print("Extracting taxid 2 taxonomy map for assemblies of interest using taxonkit")
        os.system("cut -f 6 " + assemblypath + " > taxids.txt")
        os.system("./taxonkit --data-dir . lineage -t taxids.txt > taxonkit_output.txt")

    # parse the taxonkit output into ncbi ID : taxonomy
    ncbi2tax = parse_taxonkit_output('taxonkit_output.txt') # taxid : taxonomy
    acc2tax = {} # accession number : taxonomy
    
    # loop through assembly file, gather ftp link and taxonomy
    print("Building assembly 2 taxonomy map")
    with open(assemblypath,'r') as f:
        for line in f:
            words = line.strip().split('\t')
            acc = words[0]
            taxid = words[5]
            acc2tax[acc] = ncbi2tax[taxid]

    # remove intermediate files
    to_remove = ['names.dmp','nodes.dmp','delnodes.dmp','merged.dmp','taxids.txt','taxdump.tar.gz']
    for ff in to_remove:
        if os.path.exists(ff):
            os.remove(ff)

    # write taxonomy map
    if outfile is not None:
        print('Writing taxonomy map')
        with open(outfile,'w') as f:
            keys = sorted(acc2tax.keys())
            for key in keys:
                f.write(key + '\t' + acc2tax[key] + '\n')
    else:
        return acc2tax

# uses taxonkit tool to download and make mapping of
# RefSeq accession to taxonomy
# 
def get_locus2taxonomy(assemblypath,fnapath,delim="|",outfile=None):
    """
    Makes a mapping of fna locus IDs to taxonomy strings.
    
    Note: The output tax file is for a gene-split database
    mapping loci (e.g. GCF_000007365.1|WP_011053539.1) to taxonomy.
    Assumes headers are delimited with accession first and 
    gene ID second.

    Parameters
    ----------
    assemblypath : str
        The file location of the assembly summary from refseq
        containing the accessions of interest in column one and
        the NCBI taxon ID in column 6.

    fnapath : str
        The fna database file containing all relevant loci.
        
    outfile : str, optional
        Output path for tab-delimited tax file.

    Returns
    -------
    dict
        A dict mapping accession:taxonomy (or nothing if
        outfile is passed).
    """

    # get accession 2 taxonomy
    acc2tax = get_accession2taxonomy(assemblypath)
    locus2tax = {} # locusID : taxonomy
    # loop through fna file, gather taxonomy map
    print("Building locus 2 taxonomy map")
    count = 0
    with open(fnapath,'r') as f:
        for line in f:
            if line.startswith('>'):
                count += 1
                header = line[1:].split()[0] # drop comments and ">"
                words = header.split(delim) # split on delimiter
                acc = words[0]
                if header in locus2tax:
                    print("Warning: header already found previously: " + header)
                try:
                    locus2tax[header] = acc2tax[acc]
                except KeyError:
                    print('Warning: accesion " + header + " not found in acc2tax.')
    print('Processed ' + str(count) + ' DNA sequences')
    print('There are ' + str(len(locus2tax)) + ' keys in the locus2tax dict.')
    
    # write taxonomy map
    if outfile is not None:
        print('Writing taxonomy map')
        with open(outfile,'w') as f:
            keys = sorted(locus2tax.keys())
            for locus in keys:
                f.write(locus + '\t' + locus2tax[locus] + '\n')
    else:
        return locus2tax

# parses taxon kit output into ncbiID : taxonomy mapping
#
# writes to output file if provided, or returns a dict
def parse_taxonkit_output(taxonkit_output,outfile=None):
    id2tax = {}
    with open(taxonkit_output,'r') as f:
        for line in f:
            cols = line.strip().split('\t')
            taxid = cols[0]
            taxraw = cols[1]
            taxraw = taxraw.split(';')
            tax = 'k__' + taxraw[1]
            tax += ';p__'
            if len(taxraw) > 2:
                tax += taxraw[2]
            tax += ';c__'
            if len(taxraw) > 3:
                tax += taxraw[3]
            tax += ';o__'
            if len(taxraw) > 4:
                tax += taxraw[4]
            tax += ';f__'
            if len(taxraw) > 5:
                tax += taxraw[5]
            tax += ';g__'
            if len(taxraw) > 6:
                tax += taxraw[6]
            tax += ';s__'
            if len(taxraw) > 7:
                tax += taxraw[7]
            tax += ';t__'
            if len(taxraw) > 9:
                if len(taxraw[9]) > 0:
                    tax += taxraw[9]
                else:
                    tax += taxraw[7]
            elif len(taxraw) == 9:
                if len(taxraw[8]) > 0:
                    tax += taxraw[8]
                else:
                    tax += taxraw[7]
            else:
                tax += 'None'
            id2tax[taxid] = tax

    if outfile is None:
        return id2tax
    else:
        with open(outfile,'w') as f:
            keys = sorted(id2tax.keys())
            for key in keys:
                f.write(key + '\t' + id2tax + '\n')

# download fasta-formatted gene-split genomes
# and taxonomy file
def make_refseq_fasta_and_taxonomy(assemblypath, dbpath, taxpath):
    donelist = set()  # list of completed accessions
    outdir = os.path.dirname(dbpath)
    
    # make dir or load partial db
    if os.path.exists(dbpath):
        print("Existing database found at " + dbpath + ". Loading...")
        with open(dbpath,'r') as f:
            for line in f:
                if line.startswith('>'):
                    acc = '_'.join(line[1:].split('_')[:2])
                    donelist.add(acc)
        print("Found " + str(len(donelist)) + " existing genomes to skip.")
    else:
        if not os.path.exists(outdir):
            print("Making output directory " + outdir)
            try:
                os.makedirs(outdir)
            except OSError:
                print ("Refusing to overwrite output dir %s" % outdir)
                raise

    # loop through assembly file, gather ftp links
    ftplinks = {} # refseq accession:ftp link
    print("Loading ftp links for assemblies")
    with open(assemblypath,'r') as f:
        for line in f:
            words = line.strip().split('\t')
            acc = words[0]
            ftp = words[19]
            ftplinks[acc] = ftp

    # for each strain, download file, rename headers with accID
    # example header: >lcl|NC_004061.1_cds_WP_011053539.1_1
    # make it pipe-delimited, e.g. >GCF_000010525.1|WP_011053539.1
    count = 0
    print('Downloading and processing genomes')
    with open(dbpath,'a') as f:
        keys = sorted(ftplinks.keys())
        for acc in keys:
            count += 1
            sys.stdout.write(str(count) + '/' + str(len(ftplinks)) + ' ')
            sys.stdout.flush()
            if acc in donelist:
                continue
            basename = os.path.basename(ftplinks[acc])
            filename = basename + '_cds_from_genomic.fna.gz'
            os.system("wget --retry-connrefused --waitretry=1 --read-timeout=20 --timeout=15 -t 0 -O " + filename + " " + ftplinks[acc] + '/' + filename + " >& /dev/null ")
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
                        comments = line.strip()[(line.index(' ')+1):] # comments after whitespace

                        # example before: >lcl|NC_004061.1_cds_WP_011053539.1_385
                        # this is >lcl|<ncbi id>_cds_<refseq geneID>_<gene index>
                        # sometimes there is no refseq geneID:
                        # >lcl|<ncbi id>_cds_<gene index>
                        # make it pipe-delimited, e.g. >GCF_000010525.1|WP_011053539.1|385
                        # or, e.g. >GCF_000010525.1||385
                        # Note: it is necessary to keep the gene index
                        # because some genes show up twice
                        # 1. get everything after "_cds_"
                        genedesc = header[(header.index('_cds_')+5):]
                        # 2. if there are "_" in the gene description,
                        #    split the gene ID from the gene index at end.
                        #    otherwise geneID is empty
                        if '_' in genedesc:
                            geneID = genedesc[:genedesc.rindex('_')]
                            geneindex = genedesc[(genedesc.rindex('_')+1):]
                        else:
                            geneID = ''
                            geneindex = genedesc
                        # 3. join the fields with delimiter "|"
                        header = '>' + '|'.join([acc,geneID,geneindex])                        
                        header += ' ' + comments # add back comments at end
                    else:
                        seq += line.strip()
            f.write(header + '\n' + seq + '\n') # don't forget to write the last sequence
            os.remove(filename)
        f.flush() # force write of current genome to output file
        sys.stdout.write('\n')

    print('Finished downloading ' + str(count) + ' genomes.')

    # create and write taxonomy file if it doesn't exist
    if not os.path.exists(taxpath):
        print('Creating taxonomy file.')
        get_locus2taxonomy(assemblypath,dbpath,outfile=taxpath)
    else:
        print("Taxonomy file " + taxpath + " already exists; skipping creation.")
