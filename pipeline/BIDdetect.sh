#!/bin/bash
  
ml R

# Get the directory where BIDdetect.sh is located
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# First argument is directory with sorted bam files NO FINAL SLASH
# Second argument is where outputs should go
# Third argument is path to reference fasta file
# Fourth argument is path to bed file with regions to examine

# Base path to the sorted & indexed bam files
BAM_PATH="$1"
DEST_PATH="$2"

# Check if BAM_PATH exists
if [ ! -d "$BAM_PATH" ]; then
    echo "Error: BAM directory '$BAM_PATH' does not exist."
    exit 1
fi

# Create DEST_PATH if it doesn't exist
if [ ! -d "$DEST_PATH" ]; then
    echo "Creating destination directory: $DEST_PATH"
    mkdir -p "$DEST_PATH"
    if [ $? -ne 0 ]; then
        echo "Error: Failed to create destination directory '$DEST_PATH'."
        exit 1
    fi
fi 

# Name other input files
ref_fa="$3"
bed="$4"

# Check if reference fasta and bed files exist
if [ ! -f "$ref_fa" ]; then
    echo "Error: Reference fasta file '$ref_fa' does not exist."
    exit 1
fi

if [ ! -f "$bed" ]; then
    echo "Error: BED file '$bed' does not exist."
    exit 1
fi

# Name where final counts will go
count_file="${DEST_PATH}/BIDdetect_counts.txt"
# Log file location
LOG_FILE="${DEST_PATH}/BIDdetect.log"

# Initialize the final output file with headers
echo -e "sample\tchr\tpos\tgene\ttotalReads\tA.count\tC.count\tG.count\tT.count\tDeletion.count\tInsertion.count\tref\tkmer\tstrand\tdelrate" > "$count_file"

# Create or clear the log file
echo "Starting BIDdetect processing..." > "$LOG_FILE"
echo "Log file created at: $LOG_FILE" >> "$LOG_FILE"

(cd "$BAM_PATH"
# Loop through each .bam file in the input directory
for bam_file in "$BAM_PATH"/*.bam; do

    sample=$(basename "${bam_file}" _sort.bam)  # Get the base name without the extension

    echo "Processing sample: ${sample}" #>> "$LOG_FILE"

    # Extract counts using bam_counts.R
    echo "Running R script to count for sample: ${sample}" >> "$LOG_FILE"
    Rscript "${SCRIPT_DIR}/bam_counts.R" --bedFile "$bed" --bamFile "$bam_file" --referenceFasta "$ref_fa" --outputFile "${sample}_counts.txt" >> "$LOG_FILE" 2>&1

    # Process the resultant file and append it to the final output
    result_file="${sample}_counts.txt"

    # Append the rest of the lines (skip header)
    awk -v sample="$sample" 'NR > 1 {print sample "\t" $0}' "$result_file" >> "$count_file"

    echo "Counts for sample ${sample} appended to ${count_file}" >> "$LOG_FILE"
done
)

echo "All files counted. Final count file located at $count_file" >> "$LOG_FILE"

output_file="${DEST_PATH}/BIDdetect_data.txt"

# Split sample name into different columns
echo "Running R script to split sample names..." >> "$LOG_FILE"
Rscript "${SCRIPT_DIR}/sample_name.R" "$count_file" "$output_file" celltype_vector_rep_treat >> "$LOG_FILE" 2>&1

echo "Processing completed. Log file located at $LOG_FILE" >> "$LOG_FILE"