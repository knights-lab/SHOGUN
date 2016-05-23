#!/usr/bin/env bash
# set up the directory structure
mkdir /project/flatiron/ben/data/ribo/trimmed_seqs

# load in the sequences
sshfs hillm096@login.msi.umn.edu:/home/knightsd/data_release/umgc/hiseq/160126_700506R_0537_BHHJT5BCXX/Knights_Project_022 /export/scratch/ben/msi

# change directories into the sequences
cd /scratch/ben/msi

# trim and clip adapters on PE data
for f in /export/scratch/ben/msi/*R1_001.fastq; do echo "bh_trimmomaticomatic $f ${f//R1/R2} 4 15 150" >>! /project/flatiron/ben/data/ribo/trimmed_seqs/trimmomaticomatic.sh; done

# run the trimmer
cd /project/flatiron/ben/data/ribo/trimmed_seqs
bash trimmomatic.sh

# make the img/silva directories
mkdir /project/flatiron/ben/data/ribo/silva_genes_aligned
mkdir /project/flatiron/ben/data/ribo/img_genes_aligned

# create the alignment files
for f in *R1_001.fastq; do echo "bowtie2 --no-unal -x /project/flatiron/data/db/fasta/bt2/SILVA_119_SSU_LSU_combined -S /project/flatiron/ben/data/ribo/silva_genes_aligned/${f//.fastq/.sam} --np 0 --mp "1,1" --rdg "0,1" --rfg "0,1" --score-min "L,0,-0.02" --norc -q -1 /project/flatiron/ben/data/ribo/trimmed_seqs/$f -2 /project/flatiron/ben/data/ribo/trimmed_seqs/${f//R1/R2} --very-sensitive -p 24 --no-hd -k 8" >>! bowtie_silva_align.sh; done
for f in *R1_001.fastq; do echo "bowtie2 --no-unal -x /project/flatiron/tonya/img_bowtie_builds/img.gene.bacteria.bowtie -S /project/flatiron/ben/data/ribo/img_genes_aligned/${f//.fastq/.sam} --np 0 --mp "1,1" --rdg "0,1" --rfg "0,1" --score-min "L,0,-0.02" --norc -q -1 /project/flatiron/ben/data/ribo/trimmed_seqs/$f -2 /project/flatiron/ben/data/ribo/trimmed_seqs/${f//R1/R2} --very-sensitive -p 24 --no-hd -k 8" >>! bowtie_img_align.sh; done

 # run the alignments
bash bowtie_silva_align.sh
bash bowtie_img_align.sh

# run lca on all sam files
cd /project/flatiron/ben/data/ribo/img_genes_aligned
python /project/flatiron/ben/NINJA-Shogun/src/bin/otu_lca.py -i $(pwd) -p 'img' -o genus.csv -d 6
python /project/flatiron/ben/NINJA-Shogun/src/bin/otu_lca.py -i $(pwd) -p 'img' -o species.csv


cd /project/flatiron/ben/data/ribo/silva_genes_aligned
python /project/flatiron/ben/NINJA-Shogun/src/bin/otu_lca.py -i $(pwd) -p 'silva' -o genus.csv -d 6
python /project/flatiron/ben/NINJA-Shogun/src/bin/otu_lca.py -i $(pwd) -p 'silva' -o species.csv

for f in *R1_001.fastq; do grep '^@' ${f} | wc -l >>! num_reads.txt; done
