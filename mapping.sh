#!/bin/bash

ml biology samtools
ml biology py-cutadapt/1.18_py36

# first argument $1 is location of fastq files NO FINAL SLASH
DIR_fastq = $1
# second argument $2 is location of template aligning to
ref_fa = $2
# third arugment $3 is output location for mapped files NO FINAL SLASH
DIR_dest = $3

(cd $1
# Loop through each .fq file in the input directory
for fq_file in "$1"/*.fq; do
    sample=$(basename "${fq_file}" .fq)  # Get the base name without the extension

    (
        echo "Processing $sample..."
        
        echo "removing adapter"
        # trim adapters with CutAdapt with all cores available
        cutadapt -a CACTCGGGCACCAAGGAC -g GACGCTCTTCCGATCT -o $DIR_fastq/trimmed/"${sample}"_trimmed.fq "${fq_file}"

        echo "mapping with minimap2"
        # Map using minimap2
        /oak/stanford/groups/nicolemm/rodell/minimap2/minimap2 -a $ref_fa.fasta $DIR_fastq/trimmed/"${sample}"_trimmed.fq -k5 -t 15 > $DIR_dest/minimap2/"${sample}.sam"

        # Convert to BAM file, sort the file, and index it
        echo "Converting and sorting for $sample mapped with minimap2"
        samtools view -b -o "$DIR_dest/minimap2/${sample}.bam" "$DIR_dest/minimap2/${sample}.sam"
        samtools sort -@ 15 -O bam -o "$DIR_dest/minimap2/sorted/${sample}_sort.bam" "$DIR_dest/minimap2/${sample}.bam"
        (cd "$DIR_dest/minimap2/sorted/" && samtools index "${sample}_sort.bam")

        if ls $DIR_dest/minimap2/sorted/"${sample}_sort.bam"; then
                echo "$sample mapping with minimap2 complete."

        else
            echo "$sample mapping with minimap2 failed. No sorted bam file exists."
        fi
    )
done )

wait # Wait for all background jobs to finish
