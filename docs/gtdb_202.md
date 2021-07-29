## Dependencies
* helix
* BURST -> as well as the embalmlet binaries
* seqkit

## GTDB release 06-R202

### Download and extract the files
```
# Download the representative genomes
wget --retry-connrefused --waitretry=1 --read-timeout=20 --timeout=15 -t 40 https://data.gtdb.ecogenomic.org/releases/release202/202.0/genomic_files_reps/gtdb_genomes_reps_r202.tar.gz

# Extract the genomes
tar -xzvf ./gtdb_genomes_reps_r202.tar.gz
```

```
# Move all files to current folder
for f in **/*.gz; do mv $f ./; done

# Strip trailing characters from the filenames
for f in *.fna.gz; do FN=${f/_/-}; FN=${FN/_*/}; mv $f ${FN/-/_}.fna.gz; done

# Linearize the database
/mnt/nvidia/pkr/code/BURST/bin/lingenome . combined_seqs.fna FILENAME
```

### Build the taxonomy files 
```
cd ../
mkdir taxonomy
cd taxonomy

# Download and sort taxonomy file
wget https://data.ace.uq.edu.au/public/gtdb/data/releases/release202/202.0/ar122_taxonomy_r202.tsv -O- >> taxonomy.tsv
wget https://data.ace.uq.edu.au/public/gtdb/data/releases/release202/202.0/bac120_taxonomy_r202.tsv -O- >> taxonomy.tsv

# Extract the headers from a fasta
grep -e "^>" ../combined_seqs.fna > headers.txt

cut -c 4- < taxonomy.tsv > taxonomy.stripped.tsv

/usr/bin/time -v python /mnt/nvidia/pkr/code/helix/helix/gtdb_taxonomy.py -f ./headers.txt -t ./taxonomy.tsv -o taxonomy.tax
```

### Build the BURST database
```
# Build the database (on MSI must have greater than 500 GB RAM)
/usr/bin/time -v /panfs/roc/groups/8/knightsd/hillm096/BURST/bin/burst_linux_DB15 -r ./combined_seqs.fna -o ./combined_seqs.s3333.d301.edx -a ./combined_seqs.s3333.d301.acx -d DNA 301 -i 0.95 -t 128 -s 3333
```

### Building the shear_bayes.txt
``` 
# Shear the database
# This is from HELIX 
/usr/bin/time -v python /mnt/nvidia/pkr/code/helix/helix/shear_db -f ./combined_seqs.fna -r 100 -s 100 -o ./shear_100_100.fna

# align to the database
/usr/bin/time -v /panfs/roc/groups/8/knightsd/hillm096/BURST/bin/burst_linux_DB15 -t 128 -q ./shear_100_100.fna -a ./combined_seqs.s3333.d301.acx -r combined_seqs.s3333.d301.edx -o ../burst_shear.b6 -m ALLPATHS -i 0.98;

sed 's/_/./1' burst_shear.b6 | sed 's/_/./1' > burst_shear.fixed.b6

/mnt/nvidia/pkr/code/BURST/embalmlets/bin/embalmulate burst_shear.fixed.b6 sheared_bayes.fixed.otu.txt

/usr/bin/time -v python /mnt/nvidia/pkr/code/helix/helix/shear_results.py -t ./sheared_bayes.fixed.otu -m ./taxonomy/taxonomy.tax --output ./sheared_bayes.txt
```

# Genome annotation

For genome annotation, we had to download and process the genes file provided by the GTDB.
Genes were identified with prodigal for the GTDB internal process.
We annotated the genes using the eggnog mapper version 2.1.4 from install from Bioconda.

```
emapper.py -m diamond --no_annot --no_file_comments -i combined_seqs.fna -o annotations/archaea --data_dir /mnt/btrfs/data/eggnog --itype proteins --tax_scope Archaea --cpu 40
emapper.py  --annotate_hits_table archaea.emapper.seed_orthologs --no_file_comments -o archaea --cpu 40 --dbmem
```

Archaea annotation did not take that long using Diamond.
```
================================================================================

Total hits processed: 511114
Total time: 4788 secs
FINISHED
```


Bacteria gene annotation took 8 days on a 96 core machine.
```
emapper.py -m diamond --no_annot --no_file_comments -i combined_seqs.fna -o bacteria --data_dir /project/flatiron2/data/eggnog --itype proteins --tax_scope Bacteria --cpu 96
emapper.py  --annotate_hits_table bacteria.emapper.seed_orthologs --no_file_comments -o bacteria --cpu 96 --dbmem --data_dir /project/flatiron2/data/eggnog
```

``````
================================================================================

Total hits processed: 135279611
Total time: 693713 secs
FINISHED

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
