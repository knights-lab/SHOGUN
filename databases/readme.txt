# step 1. download refseq "GMG" gene fasta and taxa using Gabe's scripts
# file will be named GMG.bacterial.genes.fasta

# step 2. scrape protein IDs and descriptions from fasta headers
# script is in the SHOGUN scripts folder
python shogun/shogun/scripts/scrape_organism2protein_map_from_GMG.py GMG.bacterial.genes.fasta GMG.bacteria.to.protein.ID.map.txt GMG.bacteria.protein.ID.to.description.map.txt

# step 3. download uniprot mapping from here
wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/idmapping.dat.gz

# step 4. parse the refseqID 2 KEGG KO ID from this idmapping file
python /project/flatiron2/dan/shogun/shogun/scripts/make_refseq_2_KO_map_from_uniprot_ID_mapping.py idmapping.dat rs2ko.txt

# step 5. convert the organism 2 refseq protein mapping to an org 2 kegg mapping
time python /project/flatiron2/dan/shogun/shogun/scripts/convert_GMG.bacteria.protein.table.from.refseq.to.x.py GMG.bacteria.to.protein.ID.map.txt rs2ko.txt GMG.bacteria.to.KO.map.txt


# optional: download eggnog hierarchy from here? not sure how to use this
wget http://eggnogdb.embl.de/download/eggnog_4.5/OG_hierarchies.tsv.gz

