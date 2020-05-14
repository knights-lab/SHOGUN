# Gene split database

Download the newest idmapping.data.gz from UniProt ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/

https://raw.githubusercontent.com/knights-lab/analysis_SHOGUN/master/scripts/shear_db.py

```bash
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/assembly_summary.txt -O- > assembly_summary.txt
awk -F "\t" '($5 == "representative genome" || $5 == "reference genome") && $14=="Full" && $11=="latest"{print $20}' assembly_summary.txt | head -n 100 > ftpdirpaths
awk 'BEGIN{FS=OFS="/";filesuffix="cds_from_genomic.fna.gz"}{ftpdir=$0;asm=$10;file=asm"_"filesuffix;print ftpdir,file}' ftpdirpaths > ftpfilepaths
cat ftpfilepaths | xargs -n 1 -P 16 wget -q --retry-connrefused --waitretry=1 --read-timeout=20 --timeout=15 -t 99 -P downloads/
gunzip ./downloads/*

# Grab the lingenome
wget https://github.com/knights-lab/BURST/raw/master/embalmlets/bin/lingenome
chmod +x ./lingenome
./lingenome downloads combined_seq.fna combined_seqs.plasmids.fna HEADFIX
```

Note that some of the genes have pseudo=True and won't have annotation in UniProt.

```python

```

```bash
wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/idmapping.dat.gz
gunzip idmapping.data.gz
```

```
mkdir burst
burst15 --references genomes.small.fna --taxonomy genomes.small.tax --output burst/genomes.small.edx --npenalize --makedb --fingerprint --accelerator burst/genomes.small.acx
```

```
# How to make sheared_bayes.txt
python shear_db.py Rep94.fasta 100 50 > Rep94.shear.fasta
/project/flatiron2/sop/burst15 -t 16 -q Rep94.shear.fasta -a Rep94.acx -r Rep94.edx -o Rep94.shear.b6 -m ALLPATHS -fr -i 0.98
sed 's/_/./1' Rep94.shear.b6 > Rep94.shear.fixed.b6
/project/flatiron2/sop/embalmulate Rep94.shear.fixed.b6 Rep94.shear.otu.txt
python shear_results_fix.py Rep94.shear.otu.txt Rep94.tax Rep94.shear
```


https://raw.githubusercontent.com/knights-lab/analysis_SHOGUN/master/scripts/shear_db.py



```
# Make database for UTree
mkdir utree
utree-build_gg genomes.small.fna genomes.small.tax genomes.small.utr 8 0 RC

```

```
# Shear the database
```

```
# Prepare for Prokka
```
