#!/bin/bash
#
# This script first aggregates all '*.run_metrics.tsv' files from a pipeline run
# into a single combined file. It then parses that file to extract the read counts
# for the deduplicated BAM step, creating a final, clean summary table.
#

# Stop the script if any command fails or if an unset variable is used
set -euo pipefail

# --- 1. Define File and Directory Paths ---
# The main directory where all the sample subdirectories are located.
TOP_LEVEL_OUTPUT_DIR="/scratch/users/rodell/20251022_endoBID/mapping"

# The intermediate file where all metrics will be combined.
COMBINED_METRICS_FILE="${TOP_LEVEL_OUTPUT_DIR}/combined_run_metrics.tsv"

# The final, two-column summary file.
DEDUP_COUNTS_FILE="${TOP_LEVEL_OUTPUT_DIR}/deduplicated_counts_summary.tsv"


# --- 2. Aggregate All Individual Metrics Files ---
echo "--- Step 1: Aggregating all run_metrics.tsv files ---"

# This command finds all metrics files, grabs the header from the first one,
# and then appends the data (without headers) from all of them.
(head -n 1 $(find "$TOP_LEVEL_OUTPUT_DIR" -name "*.run_metrics.tsv" -print -quit) && \
 find "$TOP_LEVEL_OUTPUT_DIR" -name "*.run_metrics.tsv" -exec tail -n +2 {} +) > "$COMBINED_METRICS_FILE"

# Validation: Check if the combined file was created and is not empty.
if [ ! -s "$COMBINED_METRICS_FILE" ]; then
    echo "Error: Aggregation failed. The combined metrics file is empty or was not created."
    echo "Please check if '*.run_metrics.tsv' files exist in '$TOP_LEVEL_OUTPUT_DIR'."
    exit 1
fi

echo "Aggregation complete. Combined file is at: $COMBINED_METRICS_FILE"
echo ""


# --- 3. Parse the Combined File to Create the Final Summary ---
echo "--- Step 2: Parsing combined file to extract deduplicated counts ---"

awk '
    # This block runs once at the beginning to print the header.
    BEGIN { 
        printf "Sample_Name\tDeduplicated_Read_Count\n" 
    } 
    # This block runs on every line where the first column ($1) is "7_Deduplicated_BAM".
    $1 == "7_Deduplicated_BAM" { 
        # It prints the second column ($2) and fourth column ($4), separated by a tab.
        printf "%s\t%s\n", $2, $4 
    }
' "$COMBINED_METRICS_FILE" > "$DEDUP_COUNTS_FILE"

echo "Parsing complete."


# --- 4. Final Confirmation ---
echo ""
echo "--- Success! ---"
echo "Final summary of deduplicated counts is at: $DEDUP_COUNTS_FILE"
echo ""
echo "Preview of the final file:"
head "$DEDUP_COUNTS_FILE"