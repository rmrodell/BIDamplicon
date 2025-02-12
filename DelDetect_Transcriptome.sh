#!/bin/bash

# first argument is file path to the sorted & indexed bam file
# second argument is where outputs should go
# third argument is ideally formatted sample name
# fourth argument is path to reference fasta file
# fifth argument is path to bed file with regions to examine

# Base path to the sorted & indexed bam files
BAM_FILE="$1"
# Destination directory where output files should go
DEST_PATH="$2"
# name of the sample
sample="$3"
# Path to the reference fasta file
ref_fa="$4"

# Path to the BED file to be counting over
bed="$5"
# extract the interval number from the bed file
interval=$(echo "$bed" | awk -F'[_|.]' '{print $(NF-1)}')

# name where final counts will go
count_file="${DEST_PATH}/${sample}_${interval}_counts.txt"

echo "Processing sample: ${sample}"
echo "Final counts will be available at: ${count_file}"

# extract counts using bam_counts.R
Rscript bam_count.R --bedFile $bed --bamFile $bam_file --referenceFasta $ref_fa --outputFile $countfile

echo "Counting complete. Final count file located at $count_file" 

# concatenate the different count files later