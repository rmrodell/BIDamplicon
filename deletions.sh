#!/bin/bash

# first argument is directory with find sorted bam files NO FINAL SLASH
# second argument is tab separated file with combos for processing
# third argument is path to reference fasta file
# fourth argument is path to bed file with regions to examine
# fifth argument is read depth cutoff
# sixth argument is mismatch percent cutoff



# Base path to the sequencing data
BASE_PATH="$1"
FINAL_OUTPUT_FILE="${BASE_PATH}/all_results.csv"

# Name of the CSV file
combos_file="$2"

# Check if the input file exists
if [ ! -f "$combos_file" ]; then
    echo "Error: File $combos_file not found!"
    exit 1
fi

(cd /scratch/users/rodell/ModDetect
# run ModDetect on files
# Initialize the final output file
echo "sample,chr,pos,NReads_test,gene,metadata,A.count,C.count,G.count,T.count,Deletion.count,Insertion.count,A.count.ctrl,C.count.ctrl,G.count.ctrl,T.count.ctrl,Deletion.count.ctrl,Insertion.count.ctrl,reference,kmer,strand,NReads_ctrl,mm.perc,ctrl.err,expected.err,p" > "$FINAL_OUTPUT_FILE"

# Variable to track if the header has been written in final output file
header_written=false

# Read the CSV file line by line
while IFS=$'\t' read -r f_file g_file sample
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
        -r "$3" \
        -m 1 \
        -o "${BASE_PATH}/${sample}.csv" \
        -b $4 \
        -x $5 \
        -y $6 \
        -d

    # Process the resultant file and append it to the final output
    result_file="${BASE_PATH}/${sample}.csv"

    if [ -f "$result_file" ]; then

    # If the header hasn't been written yet, we write the header line
        if [ "$header_written" = false ]; then
            # Print the header of the first file and set the flag
            awk -v sample="$sample" 'NR==1 {print $0}' "$result_file" >> "$FINAL_OUTPUT_FILE"
            header_written=true
        fi
        
        # Append the rest of the lines (skip header for subsequent files)
        awk -v sample="$sample" 'NR > 1 {print sample "," $0}' "$result_file" >> "$FINAL_OUTPUT_FILE"

    else
        echo "Warning: Expected output file $result_file not found!"
    fi

done < "$combos_file")

echo "All files processed with ModDetect. Final summary file located at ${BASE_PATH}/all_results.csv" 