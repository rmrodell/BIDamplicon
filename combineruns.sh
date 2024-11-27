#!/bin/bash

# Define the source directories
dir1="/scratch/users/rodell/20241114_pool1/BR_pool1/20241112/fastq"
dir2="/scratch/users/rodell/20241114_pool1/BR_pool1/20241113/fastq"
target_dir="/scratch/users/rodell/20241114_pool1/BR_pool1/fastq"

# Loop through each FASTQ file in the first directory
for file in "$dir1"/*.fq; do
    # Extract the base filename (without path)
    filename=$(basename "$file")
    echo "Processing $filename"
    # Check if the corresponding file exists in the second directory
    if [[ -f "$dir2/$filename" ]]; then
        # Concatenate files and write to the target directory
        cat "$file" "$dir2/$filename" > "$target_dir/$filename"
        echo "Concatenated $filename -> $target_dir/$filename"
    else
        echo "Warning: $dir2/$filename not found"
    fi
done
