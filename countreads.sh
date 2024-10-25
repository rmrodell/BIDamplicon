#!/bin/bash

ml biology samtools

# Directory containing the SAM files
SAM_DIR="/oak/stanford/groups/nicolemm/rodell/BIDamplicon/minimap2/20241003"
# Output file for the count results
output_file="read_counts.txt"

# Create or clear the output file
echo "SAM File Name, Read Count" > "$output_file"

# Iterate over all .sam files in the specified directory
for sam_file in "$SAM_DIR"/*.sam; do
    # Check if the file exists
    if [[ -f "$sam_file" ]]; then
        # Count the number of reads using samtools
        read_count=$(samtools view -c -F 4 "$sam_file")

        # Get the base name of the SAM file
        base_name=$(basename "$sam_file")

        # Append the result to the output file
        echo "$base_name, $read_count" >> "$output_file"
    else
        echo "No .sam files found in $SAM_DIR."
    fi
done

echo "Read counts saved to $output_file."