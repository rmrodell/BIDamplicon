#!/bin/bash

ml biology samtools
ml R

# first argument $1 should be rest of file directory to source location of fastq files NO FINAL SLASH
# second argument $2 should be rest of file directory to final location of fq files NO FINAL SLASH
# third argument $3 is the file.csv containing barcodes and filenames
# fourth argument $4 is the directory to the file you are aligning your sequences to
# fifth argument $5 is directory to the file with the combos ModDetect will be run on
# sixth argument $6 is the bed file ModDetect will use to find sites in

# Define the base directory containing your fastq files from the Nanopore
BASE_DIR="/oak/stanford/groups/nicolemm/rodell/BIDamplicon/$1"
# Define the destination directory for the concatenated files
DEST_DIR="/oak/stanford/groups/nicolemm/rodell/BIDamplicon/$2"

# concatenate files, rename to the actual file name, move to be within the minimap2 destination directory

# Set the input CSV file
barcode_file="$3"

# Read the CSV file line by line
while IFS=, read -r barcode filename
do
    # Skip empty lines
    [ -z "$barcode" ] || [ -z "$filename" ] || {
        # Call your original command here; this is just an example
        echo "Processing barcode: $barcode with filename: $filename"
        
        # Navigate to the directory corresponding to the barcode
    DIR="$BASE_DIR/$barcode"
    
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
        cp "$BASE_DIR/$barcode/$output_fastq" "$DEST_DIR/"

        # Generate MD5 checksum for the concatenated file in the destination directory
        md5_dest=$(md5sum "${DEST_DIR}/${output_fastq}" | awk '{ print $1 }')
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

(cd $DEST_DIR

# align sequences

# Read the CSV file line by line
while IFS=, read -r barcode filename
do
    # Skip empty lines
    [ -z "$barcode" ] || [ -z "$filename" ] || {
        # Align to reference sequence using minimap2
        echo "Aligning $filename..."
        /oak/stanford/groups/nicolemm/rodell/minimap2/minimap2 -a $4.fas "$filename".fq -k5 > "$filename".sam
        
        # Convert to BAM file, sort file, index file
        echo "Converting and sorting for $filename..."
        samtools view -b -o "$filename".bam "$filename".sam
        samtools sort -O bam -o "sorted/${filename}_sort.bam" "$filename".bam
        cd "sorted"
        samtools index "${filename}_sort.bam"
        
        echo "$filename mapping complete."
    }
done < "$barcode_file")

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

# run ModDetect on files
(cd  /oak/stanford/groups/nicolemm/rodell/BIDamplicon/ModDetect/

# Base path to the sequencing data
BASE_PATH="$DEST_DIR/sorted/"
FINAL_OUTPUT_FILE="${BASE_PATH}/all_results.csv"

# Name of the CSV file
combos_file="$5"

# Check if the input file exists
if [ ! -f "$combos_file" ]; then
    echo "Error: File $combos_file not found!"
    exit 1
fi

# Initialize the final output file
echo "sample,chr,pos,NReads_test,gene,metadata,A.count,C.count,G.count,T.count,Deletion.count,Insertion.count,A.count.ctrl,C.count.ctrl,G.count.ctrl,T.count.ctrl,Deletion.count.ctrl,Insertion.count.ctrl,reference,kmer,strand,NReads_ctrl,mm.perc,ctrl.err,expected.err,p" > "$FINAL_OUTPUT_FILE"

# Read the CSV file line by line
while IFS=, read -r f_file g_file sample
do
    # Skip header or empty lines
    if [[ "$f_file" == "f_file" || -z "$f_file" ]]; then
        continue
    fi

    echo "Processing files: treated $f_file, control $g_file, sample $sample"

    # Rscript command with values from the CSV
    Rscript ModDetect.R \
        -f "${BASE_PATH}/${f_file}" \
        -g "${BASE_PATH}/${g_file}" \
        -r "$4" \
        -m 1 \
        -o "${BASE_PATH}/${sample}.csv" \
        -b "${BASE_PATH}/$6" \
        -x 1 \
        -y 1 \
        -a

    # Process the resultant file and append it to the final output
    result_file="${BASE_PATH}/${sample}.csv"

    if [ -f "$result_file" ]; then
        # Skip the header in subsequent files and prepend the sample name
        # Use awk to add the sample name to each line after the header
        awk -v sample="$sample" 'NR==1 {print $0} NR>1 {print sample "," $0}' "$result_file" >> "$FINAL_OUTPUT_FILE"
    else
        echo "Warning: Expected output file $result_file not found!"
    fi

done < "$combos_file")

echo "All files processed with ModDetect. Final summary file located at ${BASE_PATH}/all_results.csv" 

# calculate delta delta deletion rate in WT to OE or KD