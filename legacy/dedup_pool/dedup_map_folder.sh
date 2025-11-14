#!/bin/bash
set -e
set -u

# Get the directory where this script is located
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
# Get the root directory of the git repo (assuming the script is in the root)
REPO_ROOT="$SCRIPT_DIR"

# This trims two sets of adapters, so it optimized to work on pool 1 setup

# Check if required arguments are provided
if [ $# -lt 3 ]; then
    echo "Usage: $0 <input_directory> <reference_fasta> <base_directory> [threads]"
    echo "Example: $0 /path/to/fastq/dir /path/to/reference.fasta /scratch/users/username/output_dir 7"
    echo "Note: If threads not specified, default value of 1 will be used"
    exit 1
fi

# Assign command line arguments to variables
INPUT_DIR=$1
REFERENCE_FASTA=$2
BASE_DIR=$3
THREADS=${4:-1}

# Verify input directory exists
if [ ! -d "${INPUT_DIR}" ]; then
    echo "Error: Input directory ${INPUT_DIR} not found"
    exit 1
fi

# Verify reference exists
if [ ! -f "${REFERENCE_FASTA}" ]; then
    echo "Error: Reference fasta ${REFERENCE_FASTA} not found"
    exit 1
fi

# Process each .fq file in the input directory
for fastq in ${INPUT_DIR}/*.fq; do
    if [ -f "$fastq" ]; then
        echo "Processing $fastq"
        bash "${REPO_ROOT}/dedup_mapping.sh" "$fastq" "$REFERENCE_FASTA" "$BASE_DIR" "$THREADS"
    fi
done

# Create a combined text summary file
echo "Creating combined text summary..."
TEXT_SUMMARY="${BASE_DIR}/logs/combined_summary.txt"
echo "Combined Processing Summary" > ${TEXT_SUMMARY}
echo "----------------------------------------" >> ${TEXT_SUMMARY}

# Combine all individual summaries
for summary in ${BASE_DIR}/logs/*_summary.txt; do
    if [ -f "$summary" ]; then
        echo "" >> ${TEXT_SUMMARY}
        cat "$summary" >> ${TEXT_SUMMARY}
        echo "----------------------------------------" >> ${TEXT_SUMMARY}
    fi
done

# Create a combined summary table
echo "Creating combined summary table..."
TABLE_SUMMARY="${BASE_DIR}/logs/combined_summary_table.tsv"

# Create header
echo -e "Sample\tInitial_Reads\tAfter_First_Trim\tAfter_Second_Trim\tMapped_Reads\tDeduplicated_Reads\tFirst_Trim_Pct\tSecond_Trim_Pct\tMapped_Pct\tDedup_Pct\tFirst_Trim_Prev_Pct\tSecond_Trim_Prev_Pct\tMapped_Prev_Pct\tDedup_Prev_Pct" > ${TABLE_SUMMARY}

# Process each summary file and extract the numbers
for summary in ${BASE_DIR}/logs/*_summary.txt; do
    if [ -f "$summary" ]; then
        # Get sample name
        sample_name=$(basename "$summary" _summary.txt)
        
        # Extract read counts and percentages using grep and awk
        initial=$(grep "Initial reads:" "$summary" | awk '{print $3}')
        first_trim=$(grep "After first trimming:" "$summary" | awk '{print $4}')
        second_trim=$(grep "After second trimming:" "$summary" | awk '{print $4}')
        mapped=$(grep "Mapped reads:" "$summary" | awk '{print $3}')
        dedup=$(grep "After deduplication:" "$summary" | awk '{print $3}')
        
        # Extract percentages relative to initial
        first_trim_pct=$(grep "First trimming:" "$summary" | grep "initial reads" | awk '{print $3}' | sed 's/%//')
        second_trim_pct=$(grep "Second trimming:" "$summary" | grep "initial reads" | awk '{print $3}' | sed 's/%//')
        mapped_pct=$(grep "Mapped:" "$summary" | grep "initial reads" | awk '{print $3}' | sed 's/%//')
        dedup_pct=$(grep "After deduplication:" "$summary" | grep "initial reads" | awk '{print $3}' | sed 's/%//')
        
        # Extract percentages relative to previous step
        first_trim_prev=$(grep "First trimming:" "$summary" | grep "previous step" | awk '{print $3}' | sed 's/%//')
        second_trim_prev=$(grep "Second trimming:" "$summary" | grep "previous step" | awk '{print $3}' | sed 's/%//')
        mapped_prev=$(grep "Mapped:" "$summary" | grep "previous step" | awk '{print $3}' | sed 's/%//')
        dedup_prev=$(grep "After deduplication:" "$summary" | grep "previous step" | awk '{print $3}' | sed 's/%//')
        
        # Combine all values into a single line
        echo -e "${sample_name}\t${initial}\t${first_trim}\t${second_trim}\t${mapped}\t${dedup}\t${first_trim_pct}\t${second_trim_pct}\t${mapped_pct}\t${dedup_pct}\t${first_trim_prev}\t${second_trim_prev}\t${mapped_prev}\t${dedup_prev}" >> ${TABLE_SUMMARY}
    fi
done

echo
