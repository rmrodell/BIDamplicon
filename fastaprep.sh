#!/bin/bash

ml biology samtools
ml R

# first argument $1 should be file path to the first date of the sample of interest NO FINAL SLASH
# second argument $2 should be file path to the second date of the sample of interest NO FINAL SLASH
# third argument $3 is the file.csv containing barcodes and filenames
# fourth argument is the date of the first sample
# fifth argument is the date of the second sample

# Define the base directory containing your fastq files from the Nanopore
BASE_DIR="/scratch/users/rodell/20241114_pool1/BR_pool1"
# Define the destination directory for the concatenated files
# DEST_DIR="/scratch/users/rodell/20241114_pool1/BR_pool1/$1/fastq/"

# concatenate files, rename to the actual file name

# Set the input CSV file
barcode_file="$3"

# Read the CSV file line by line
while IFS=, read -r barcode filename
do
    # Skip empty lines
    [ -z "$barcode" ] || [ -z "$filename" ] || {
        # Call your original command here; this is just an example
        echo "Processing barcode: $4 $barcode with filename: $filename"
        
        # Navigate to the directory corresponding to the barcode
    DIR="$BASE_DIR/$1/$barcode"
    
    # Check if the directory exists
    if [ -d "$DIR" ]; then

        # Navigate to the directory
        cd "$DIR"

        # Unzip all .gz files
        gunzip *.gz

        echo "Files for $barcode $filename unzipped"

        # Determine the output filename for this barcode
        output_fastq="$filename.fq"

        # Concatenate all .fastq files to the specified output filename
        cat *.fastq >> "$output_fastq"

        # Move the concatenated file to the destination directory
        cp "$DIR/$output_fastq" "$BASE_DIR/$4/fastq/"

        # Generate MD5 checksum for the concatenated file in the destination directory
        md5_dest=$(md5sum "${BASE_DIR}/$4/fastq/${output_fastq}" | awk '{ print $1 }')
        md5_source=$(md5sum "${DIR}/${output_fastq}" | awk '{ print $1 }')

        # Compare the MD5 checksums
        if [ "$md5_dest" == "$md5_source" ]; then
            echo "MD5 checksum for $output_fastq matches between source and destination."
        else
            echo "WARNING: MD5 checksum for $output_fastq does NOT match between source and destination!"
        fi

        echo "Finished processing $barcode"
    else
        echo "Directory $DIR does not exist, skipping..."
    fi
    }
done 


# Read the CSV file line by line
while IFS=, read -r barcode filename
do
    # Skip empty lines
    [ -z "$barcode" ] || [ -z "$filename" ] || {
        # Call your original command here; this is just an example
        echo "Processing barcode: $5 $barcode with filename: $filename"
        
        # Navigate to the directory corresponding to the barcode
    DIR="$BASE_DIR/$2/$barcode"
    
    # Check if the directory exists
    if [ -d "$DIR" ]; then

        # Navigate to the directory
        cd "$DIR"

        # Unzip all .gz files
        gunzip *.gz

        echo "Files for $barcode $filename unzipped"

        # Determine the output filename for this barcode
        output_fastq="$filename.fq"

        # Concatenate all .fastq files to the specified output filename
        cat *.fastq >> "$output_fastq"

        # Move the concatenated file to the destination directory
        cp "$BASE_DIR/$2/$barcode/$output_fastq" "$BASE_DIR/$5/fastq/"

        # Generate MD5 checksum for the concatenated file in the destination directory
        md5_dest=$(md5sum "${BASE_DIR}/$5/fastq/${output_fastq}" | awk '{ print $1 }')
        md5_source=$(md5sum "${DIR}/${output_fastq}" | awk '{ print $1 }')

        # Compare the MD5 checksums
        if [ "$md5_dest" == "$md5_source" ]; then
            echo "MD5 checksum for $output_fastq matches between source and destination."
        else
            echo "WARNING: MD5 checksum for $output_fastq does NOT match between source and destination!"
        fi

        echo "Finished processing $barcode"
    else
        echo "Directory $DIR does not exist, skipping..."
    fi
    }
done < "$barcode_file"

echo "All barcodes have been concatenated and moved"

# concatenate files from the same barcode in a new folder

# Loop through barcode numbers from 01 to 66
for i in $(seq -f "%02g" 1 66); do
    # Construct the file names for each barcode
    file1="$BASE_DIR/$4/fastq/barcode$i.fastq"
    file2="$BASE_DIR/$5/fastq/barcode$i.fastq"
    
    # Output file name
    output_file="$BASE_DIR/fastq/barcode$i.fastq"

    # Concatenate the files if they exist
    if [[ -f "$file1" && -f "$file2" ]]; then
        cat "$file1" "$file2" > "$output_file"
        echo "Combined $file1 and $file2 into $output_file"
    else
        echo "Warning: One of the files does not exist: $file1 or $file2"
    fi
done