#!/bin/bash

ml python/3.6.1
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

echo "Processing $sample..."
        
# Extract UMIs to read name
umi_tools extract -I $file_path -p NNNNNNNNNN --log=$sample.log -S $fq_directory/{$sample}_umi.fq 

# Map using minimap2
/oak/stanford/groups/nicolemm/rodell/minimap2/minimap2 -a $ref_fa $fq_directory/{$sample}_umi.fq -k5 > $dest_dir/$sample.sam

# Convert to BAM file, sort the file, and index it
echo "Converting and sorting for $sample..."
samtools view -b -o "$dest_dir/$sample.bam" "$dest_dir/$sample.sam"
samtools sort -O bam -o "$dest_dir/sorted/{$sample}_sort.bam" "$dest_dir/${sample}.bam")
(cd "$dest_dir/sorted/" && samtools index "${sample}_sort.bam")

if ls $dest_dir/sorted/"${sample}_sort.bam"; then
    echo "$sample mapping complete."
    # deduplicate reads based on UMI
    (cd "$dest_dir/dedup/"
    umi_tools dedup -I $dest_dir/sorted/"${sample}_sort.bam" -S $dest_dir/dedup/"${sample}_dedup.bam" --log=$sample.log
    # index deduplicated file 
    samtools index "${sample}_dedup.bam"
    )
else
    echo "$sample mapping failed. No sorted bam file exists."
fi

if ls $dest_dir/dedup/"${sample}_dedup.bam"; then
    echo "$sample deduplication complete."
else
    echo "$sample deduplication failed. No deduplicated bam file exists."
fi