#!/bin/bash

ml biology samtools

# first argument $1 is directory with sorted sam files NO FINAL SLASH
# Output file for the count results
output_file="read_counts.txt"

# Create or clear the output file
echo "SAM File Name, Read Count" > "$1/$output_file"

# Iterate over all .sam files in the specified directory
for sam_file in "$1"/*.sam; do
    # Check if the file exists
    if [[ -f "$sam_file" ]]; then
        # Count the number of reads using samtools
        read_count=$(samtools view -c -F 4 "$sam_file")

        # Get the base name of the SAM file
        base_name=$(basename "$sam_file")

        # Append the result to the output file
        echo "$base_name, $read_count" >> "$1/$output_file"
    else
        echo "No .sam files found in $1."
    fi
done

echo "Read counts saved to $1/$output_file."