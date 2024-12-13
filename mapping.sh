#!/bin/bash

ml biology samtools

# first argument $1 is location of fastq files NO FINAL SLASH
# second argument $2 is location of template aligning to
# third arugment $3 is output location for mapped files NO FINAL SLASH

(cd $1
# Loop through each .fq file in the input directory
for fq_file in "$1"/*.fq; do
    filename=$(basename "${fq_file}" .fq)  # Get the base name without the extension

    (
        echo "Processing $filename..."
        
        # Map using minimap2
        /oak/stanford/groups/nicolemm/rodell/minimap2/minimap2 -a $2 "${fq_file}" -k5 -t 9 > $3/"${filename}.sam"

        # Convert to BAM file, sort the file, and index it
        echo "Converting and sorting for $filename..."
        samtools view -b -o "$3/${filename}.bam" "$3/${filename}.sam"
        samtools sort -O bam -o "$3/sorted/${filename}_sort.bam" "$3/${filename}.bam"
        (cd "$3/sorted/" && samtools index "${filename}_sort.bam")

	if ls $dest_dir/sorted/"${sample}_sort.bam"; then
    		echo "$filename mapping complete."

	else
		echo "$filename mapping failed. No sorted bam file exists."
    fi
    )
done )

wait # Wait for all background jobs to finish
