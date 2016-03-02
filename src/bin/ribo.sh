#!/usr/bin/env bash

# Set up the experiment parameters
DATA_DIR=/project/flatiron/ben/data/ribo

# Location of the sshfs from msi
INPUT_FASTQ_DIR=/export/scratch/ben/msi

#location of the fastq files
DNA_FILE_LIST=${DATA_DIR}/fastq_dna.txt
RNA_FILE_LIST=${DATA_DIR}/fastq_rna.txt

trimm_align () {
    output_dir=$1
    input_fastq_dir=$2
    file_list=$3
    bt2_indx=$4

    test -d ${output_dir} | mkdir -p  ${output_dir}
    output_fastq_dir=${output_dir}/fastq
    test -d ${output_fastq_dir} | mkdir -p  ${output_fastq_dir}
    while read in_file; do
        filename=$(basename ${in_file})
        out_file=${filename%.fastq}.trimmed.fastq
        java -jar /project/flatiron/ben/bin/Trimmomatic-0.35/trimmomatic-0.35.jar \
            SE -phred33 ${input_fastq_dir}/${in_file} \
            ${output_fastq_dir}/${out_file} \
            -threads 16 \
            LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
#            ILLUMINACLIP:/project/flatiron/ben/bin/Trimmomatic-0.35/adapters/TruSeq2-PE.fa:2:30:10 \
    done <${file_list}

    # test for directory
    output_sam_dir=${output_dir}/sam
    test -d ${output_sam_dir} | mkdir -p  ${output_sam_dir}
    # Run bowtie2
    for in_file in ${output_fastq_dir}/*.fastq; do
        # get the base filename
        filename=$(basename ${in_file})
        out_file=${filename%.fastq}.sam
        bowtie2 --no-unal --no-head -x ${bt2_indx} \
            -S ${output_sam_dir}/${out_file} --np 0 --mp "1,1" --rdg "0,1" --rfg "0,1" \
            --score-min "L,0,-0.02" --norc -q ${in_file} -k 8 -p 48 --very-sensitive
    done
}

# Copy the index into ramdisk
mkdir /dev/shm/bt2_indx/
cp /project/flatiron/gabe/IMGENES.* /dev/shm/bt2_indx
DNA_BT2_INDX=/dev/shm/bt2_indx/IMGENES
trimm_align ${DATA_DIR}/dna ${INPUT_FASTQ_DIR} ${DNA_FILE_LIST} ${DNA_BT2_INDX}
# Clear the ramdisk
rm -r /dev/shm/bt2_indx

#Copy the index into ramdisk
mkdir /dev/shm/bt2_indx/
cp /project/flatiron/data/db/fasta/bt2/SILVA_119_SSU_LSU_combined.* /dev/shm/bt2_indx
RNA_BT2_INDX=/dev/shm/bt2_indx/SILVA_119_SSU_LSU_combined
trimm_align ${DATA_DIR}/rna ${INPUT_FASTQ_DIR} ${RNA_FILE_LIST} ${RNA_BT2_INDX}
# Clear the ramdisk
rm -r /dev/shm/bt2_indx
