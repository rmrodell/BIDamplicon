#!/bin/bash

ml biology samtools

# first argument $1 should be rest of file directory to final location of fq files NO FINAL SLASH: 20241105_pool1test/analysis
# second argument $2 is the file.csv containing barcodes and filenames: /oak/stanford/groups/nicolemm/rodell/BIDamplicon/minimap2/20241105/20241105barcode.csv
# third argument $3 is the directory to the file you are aligning your sequences to:/oak/stanford/groups/nicolemm/rodell/BIDamplicon/4U.fas

# Define the directory with the concatenated fastq files
DEST_DIR="/oak/stanford/groups/nicolemm/rodell/BIDamplicon/$1"


# Set the input CSV file
barcode_file="$2"

# align sequences

# Read the CSV file line by line
while IFS=, read -r barcode filename
do
    # Skip empty lines
    [ -z "$barcode" ] || [ -z "$filename" ] || {
        # Align to reference sequence using minimap2
        echo "Aligning $filename..."
        /oak/stanford/groups/nicolemm/rodell/minimap2/minimap2 -a "$3" "$filename".fq -k5 > "$filename".sam
        
        # Convert to BAM file, sort file, index file
        echo "Converting and sorting for $filename..."
        samtools view -b -o "$filename".bam "$filename".sam
        samtools sort -O bam -o "sorted/${filename}_sort.bam" "$filename".bam
        (cd "sorted"
        echo "Current directory: $PWD"
        samtools index "${filename}_sort.bam")
        
        echo "$filename mapping complete."
        echo "Current directory: $PWD"
    }
done < "$barcode_file"

# count reads in alignment and mapping rate
(cd $DEST_DIR/sorted/

# Output file for the count results
output_counts_file="read_counts.txt"

# Create or clear the output file
echo "File Name, Read Count" > "$output_counts_file"

echo "Now counting reads:"

# Iterate over all .sam files in the specified directory
for bam_file in *.bam; do
    # Check if the file exists
    if [[ -f "$bam_file" ]]; then
        # Count the number of reads using samtools
        read_count=$(samtools view -c -F 4 "$bam_file")

        # Get the base name of the SAM file
        base_name="${bam_file%_sort.bam}"

        echo "Name: $base_name, Read Count: $read_count"

        # Append the result to the output file
        echo "$base_name, $read_count" >> "$output_counts_file"
    else
        echo "No .bam files found in $DEST_DIR/sorted/"
    fi
done

echo "Read counts saved to $output_counts_file.")