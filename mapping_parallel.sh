#!/bin/bash

ml biology samtools

# first argument is $file from the slurm iteration: /scratch/users/rodell/20241114_pool1/BR_pool1/fastq/HepG2_WT_3_input.fq
# second argument $ref_fa is location of template aligning to
# third arugment $dest_dir is output location for mapped files NO FINAL SLASH

file_path=$1
fq_directory=$(dirname "$file_path") # get the dest_dir the fastq file is located in
sample=$(basename "$file_path" .fq) # get the base name of the sample, previously sample

ref_fa=$2

dest_dir=$3

(cd $fq_directory
# Loop through each .fq file in the input dest_dir

echo "Processing $sample..."
        
# Map using minimap2
/oak/stanford/groups/nicolemm/rodell/minimap2/minimap2 -a $ref_fa $file_path -k5 -t 31 > $dest_dir/$sample.sam

# Convert to BAM file, sort the file, and index it
echo "Converting and sorting for $sample..."
samtools view -b -o "$dest_dir/$sample.bam" "$dest_dir/$sample.sam"
samtools sort -O bam -o "$dest_dir/sorted/{$sample}_sort.bam" "$dest_dir/${sample}.bam"
(cd "$dest_dir/sorted/" && samtools index "${sample}_sort.bam")
)

if ls $dest_dir/sorted/"${sample}_sort.bam"; then
    echo "$sample mapping complete."
else
    echo "$sample mapping failed. No sorted bam file exists"
fi
