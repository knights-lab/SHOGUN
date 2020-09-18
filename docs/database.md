## Dependencies
* helix
* BURST -> as well as the embalmlet binaries
* seqkit

## RefSeq release 201

### Download the files
```
# Archaea/Bacteria Database Building
# Download the assembly summary for each of the files
wget https://ftp.ncbi.nlm.nih.gov/genomes/refseq/archaea/assembly_summary.txt -O- > assembly_summary.txt
wget https://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/assembly_summary.txt -O- >> assembly_summary.txt

# Grab the latest and full genomes
awk -F "\t" '($5 == "representative genome" || $5 == "reference genome") && $14=="Full" && $11=="latest"{print $20}' assembly_summary.txt >> ftpdirpaths

# Get the ftpfilepaths from the dir paths
awk 'BEGIN{FS=OFS="/";filesuffix="genomic.fna.gz"}{ftpdir=$0;asm=$10;file=asm"_"filesuffix;print ftpdir,file}' ftpdirpaths | sed 's/ftp:/https:/' > ftpfilepaths

# Download the files
cat ftpfilepaths | xargs -n 1 -P 16 wget -q --retry-connrefused --waitretry=1 --read-timeout=20 --timeout=15 -t 40
```

### Build the taxonomy files 
```
# Strip trailing characters from the filenames
for f in *.fna.gz; do FN=${f/_/-}; FN=${FN/_*/}; mv $f ${FN/-/_}.fna.gz; done

# Linearize the database
/mnt/nvidia/pkr/code/BURST/bin/lingenome . refseq.fna FILENAME

# Download and sort taxonomy file
mkdir taxtmp && cd taxtmp
wget https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdmp.zip
unzip taxdmp.zip
t2gg nodes.dmp names.dmp tid2gg.txt SUBONLY
sort -k1,1 tid2gg.txt > tid2gg.srt.txt

# TODO comment what we are doing
cd ../fna
for f in *.fna.gz; do echo "${f/.fna.gz/}"; done | grep -F -f - ../assembly_summary.txt | cut -f 1,6,7,8 | sort -k2,2 > ../rawtax.tsv

cd ..
# TODO comment what we are doing
join -t $'\t' -12 -21 -e0 -o'1.1,2.2,1.4,0,1.3' ./rawtax.tsv taxtmp/tid2gg.srt.txt | sort -k2 > alltax.tsv

# Remove the dangling taxonomy identifiers
/usr/bin/time -v python /mnt/nvidia/pkr/code/helix/helix/strip_taxamap.py --input ./alltax.tsv --output ./taxmap.tsv
```

### Build the BURST database
```
# Build the database (on MSI must have greater than 500 GB RAM)
/usr/bin/time -v burst_linux_DB15 -r /scratch.global/ben/rep201/ab/refseq.fna -o /scratch.global/ben/rep201/ab/refseq.edx -a /scratch.global/ben/rep201/ab/refseq.acx -d DNA 301 -i 0.95 -t 32 -s 3333
```

### Building the shear_bayes.txt
``` 
# Shear the database
# This is from HELIX
# TODO 100 100 shear size
/usr/bin/time -v python shear_db -f ./gtdb.fna -r 300 -s 200 -o ./shear_300_200.fna

# align to the database
/usr/bin/time -v /panfs/roc/groups/8/knightsd/hillm096/BURST/bin/burst_linux_DB15 -t 32 -q ./shear_300_200.fna -a /scratch.global/ben/rep201/ab/refseq.acx -r /scratch.global/ben/rep201/ab/refseq.edx -o ../burst_shear.b6 -m ALLPATHS -i 0.98;

sed 's/_/./1' shear_300_200.b6 | sed 's/_/./1' > burst_shear.fixed.b6
/mnt/nvidia/pkr/code/BURST/embalmlets/bin/embalmulate burst_shear.fixed.b6 sheared_bayes.fixed.otu.txt

/usr/bin/time -v python /mnt/nvidia/pkr/code/helix/helix/shear_results.py --alignment ./sheared_bayes.fixed.otu.txt --taxa_table ../taxmap.tsv --output ./sheared_bayes.txt
```

### Creating the database.yml
```yaml
general:
  taxonomy: taxmap.tsv
  fasta: genomes.small.fna
  shear: sheared_bayes.txt
burst: burst/genomes.small
```

### Testing the database
```
# align the database
shogun align -a burst -i /project/flatiron2/ben/projects/type_1/data/hmp_mock_community/shi7_learning/combined_seqs.fna -d ./rep201_ab

# assign a rank-flexible taxonomy
shogun assign_taxonomy --no-capitalist -a burst -i ./results-200918/alignment.burst.b6 -d ./rep201_ab

# assign a rank-specific taxonomy
shogun redistribute --input ./results-200918/taxatable-200918.txt --level strain -o ./results-200918/taxatable-200918.strain.txt -d ./rep201_ab
```

## Location of the databases
```
# Location on MSI
/scratch.global/ben/shogun/rep201_ab/

# Location on teraminx
/project/flatiron/data/shogun/rep201_ab/
```


## GTDB 95
TODO
