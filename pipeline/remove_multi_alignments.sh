#!/bin/bash

# ==============================================================================
# SCRIPT: filter_primary_alignments.sh
#
# DESCRIPTION:
#   This script iterates through all BAM files in a specified input directory,
#   filters them to retain only primary alignments (removing secondary and
#   supplementary alignments), and saves the resulting clean BAM files to an
#   output directory.
#
#   It then indexes each new BAM file and generates a tab-separated summary
#   file containing the final read counts for each processed file.
#
# USAGE:
#   bash filter_primary_alignments.sh -i <input_directory> -o <output_directory>
#
# ARGUMENTS:
#   -i <input_directory>    Path to the folder containing input .bam files.
#   -o <output_directory>   Path to the folder where output files will be saved.
#
# OUTPUTS:
#   - Filtered and indexed .bam files in the output directory.
#   - A log file named 'filter_primary.log' in the output directory.
#   - A summary file named 'primary_read_counts.tsv' in the output directory.
#
# ==============================================================================

# --- Script Configuration ---
set -e
set -o pipefail
THREADS=${SLURM_CPUS_PER_TASK:-$(nproc --all)}

# --- Argument Parsing ---
while getopts "i:o:" opt; do
  case ${opt} in
    i) INPUT_DIR="$OPTARG" ;;
    o) OUTPUT_DIR="$OPTARG" ;;
    *) echo "Usage: $0 -i <input_dir> -o <output_dir>" >&2; exit 1 ;;
  esac
done

# --- Validate Arguments ---
if [ -z "${INPUT_DIR}" ] || [ -z "${OUTPUT_DIR}" ]; then
    echo "Error: Both -i (input directory) and -o (output directory) arguments are required."
    echo "Usage: bash filter_primary_alignments.sh -i <input_directory> -o <output_directory>"
    exit 1
fi

if [ ! -d "${INPUT_DIR}" ]; then
    echo "Error: Input directory '${INPUT_DIR}' not found."
    exit 1
fi

# --- Setup Output Directory and Logging ---
mkdir -p "${OUTPUT_DIR}"
LOG_FILE="${OUTPUT_DIR}/filter_primary.log"
SUMMARY_TSV="${OUTPUT_DIR}/primary_read_counts.tsv"

# Tee output to both console and log file
exec > >(tee -a "${LOG_FILE}") 2>&1

echo "=========================================================="
echo "Starting Primary Alignment Filtering"
echo "=========================================================="
echo "Input Directory:  ${INPUT_DIR}"
echo "Output Directory: ${OUTPUT_DIR}"
echo "Threads:          ${THREADS}"
echo "Log File:         ${LOG_FILE}"
echo "----------------------------------------------------------"

# --- Initialize Summary File ---
echo -e "filename\tprimary_alignment_count" > "${SUMMARY_TSV}"

# --- Main Processing Loop ---
for INPUT_BAM in "${INPUT_DIR}"/*.bam; do
    # Check if files exist to avoid errors with empty directories
    [ -f "$INPUT_BAM" ] || continue

    BASENAME=$(basename "${INPUT_BAM}")
    OUTPUT_BAM="${OUTPUT_DIR}/${BASENAME}"

    echo "Processing ${BASENAME}..."

    # Filter for primary alignments and save to new BAM
    # The flag -F 2304 excludes reads with either the secondary (256) or supplementary (2048) bit set.
    samtools view -h -b -F 2304 -@ "${THREADS}" -o "${OUTPUT_BAM}" "${INPUT_BAM}"

    # Index the newly created BAM file
    samtools index -@ "${THREADS}" "${OUTPUT_BAM}"

    # Count reads in the new file and append to summary
    COUNT=$(samtools view -c -@ "${THREADS}" "${OUTPUT_BAM}")
    echo -e "${BASENAME}\t${COUNT}" >> "${SUMMARY_TSV}"
    echo "  -> Found ${COUNT} primary alignments."

done

echo "----------------------------------------------------------"
echo "Processing complete."
echo "Summary of counts written to ${SUMMARY_TSV}"
echo "=========================================================="