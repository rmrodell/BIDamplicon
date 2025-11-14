#!/bin/bash

# Base path to the sequencing data
BASE_PATH="/oak/stanford/groups/nicolemm/rodell/BIDamplicon/minimap2/20241003"

# Arrays of specific file names
f_files=("500nMC17_WT_P102_BS_0906_sort.bam" "50nMC17_WT_P102_BS_0906_sort.bam" "50uMC1&_WT_P102_BS_0906_sort.bam" "5uMC17_WT_P102_BS_0906_sort.bam" "40nMC17_WT_P4_BS_0902_sort.bam" "4nMC17_WT_P4_BS_0902_sort.bam" "40nMC17_WT_P102_BS_0902_sort.bam" "4nMC17_WT_P102_BS_0902_sort.bam" "WT_P102_BS_0902_sort.bam" "WT_P102_BS_0906_sort.bam" "WT_P4_BS_0902_sort.bam" "WT_P4_BS_0906_sort.bam" "KD_P101_BS_0902_sort.bam" "KD_P101_BS_0906_sort.bam" "OE_P3_BS_0902_sort.bam" "OE_P3_BS_0906_sort.bam" "WT1_pLKO_BS_1214_sort.bam" "KD1_shPUS7_BS_1214_sort.bam" "WT2_pLKO_BS_1214_sort.bam" "KD2_shPUS7_BS_1214_sort.bam" "WT3_pLKO_BS_1214_sort.bam" "KD3_shPUS7_BS_1214_sort.bam" )
g_files=("500nMC17_WT_P102_in_0906_sort.bam" "50nMC17_WT_P102_in_0906_sort.bam" "50uMC1&_WT_P102_in_0906_sort.bam" "5uMC17_WT_P102_in_0906_sort.bam" "40nMC17_WT_P4_in_0902_sort.bam" "4nMC17_WT_P4_in_0902_sort.bam" "40nMC17_WT_P102_in_0902_sort.bam" "4nMC17_WT_P102_in_0902_sort.bam" "WT_P102_in_0902_sort.bam" "WT_P102_in_0906_sort.bam" "WT_P4_in_0902_sort.bam" "WT_P4_in_0906_sort.bam" "KD_P101_in_0902_sort.bam" "KD_P101_in_0906_sort.bam" "OE_P3_in_0902_sort.bam" "OE_P3_in_0906_sort.bam" "WT1_pLKO_in_1214_sort.bam" "KD1_shPUS7_in_1214_sort.bam" "WT2_pLKO_in_1214_sort.bam" "KD2_shPUS7_in_1214_sort.bam" "WT3_pLKO_in_1214_sort.bam" "KD3_shPUS7_in_1214_sort.bam")
output_files=("500nMC17_WT_P102_0906.csv" "50nMC17_WT_P102_0906.csv" "50uMC17_WT_P102_0906.csv" "5uMC17_WT_P102_0906.csv" "40nMC17_WT_P4_0902.csv" "4nMC17_WT_P4_0902.csv" "40nMC17_WT_P102_0902.csv" "4nMC17_WT_P102_0902.csv" "WT_P102_0902.csv" "WT_P102_0906.csv" "WT_P4_0902.csv" "WT_P4_0906.csv" "KD_P101_0902.csv" "KD_P101_0906.csv" "OE_P3_0902.csv" "OE_P3_0906.csv" "WT1_pLKO_1214.csv" "KD1_shPUS7_1214.csv" "WT2_pLKO_1214.csv" "KD2_shPUS7_1214.csv" "WT3_pLKO_1214.csv" "KD3_shPUS7_1214.csv" ​)
samples=("500nMC17_WT_P102_0906" "50nMC17_WT_P102_0906" "50uMC17_WT_P102_0906" "5uMC17_WT_P102_0906" "40nMC17_WT_P4_0902" "4nMC17_WT_P4_0902" "40nMC17_WT_P102_0902" "4nMC17_WT_P102_0902" "WT_P102_0902" "WT_P102_0906" "WT_P4_0902" "WT_P4_0906" "KD_P101_0902" "KD_P101_0906" "OE_P3_0902" "OE_P3_0906" "WT1_pLKO_1214" "KD1_shPUS7_1214" "WT2_pLKO_1214" "KD2_shPUS7_1214" "WT3_pLKO_1214" "KD3_shPUS7_1214" ​)  # Array of sample identifiers

# Concatenation result file
concat_output="${BASE_PATH}/output_ModDetect/20241003_psilocations.csv"

# Create/Empty the concat_output file before starting
echo "Sample,chr,pos,NReads_test,gene,metadata,A.count,C.count,G.count,T.count,Deletion.count,Insertion.count,A.count.ctrl,C.count.ctrl,G.count.ctrl,T.count.ctrl,Deletion.count.ctrl,Insertion.count.ctrl,reference,kmer,strand,NReads_ctrl,mm.perc,ctrl.err,expected.err,p" > "$concat_output"  # Adjust column headers as needed

# Loop through the arrays
for i in "${!f_files[@]}"; do
    echo "Running command for set $((i + 1))..."

    # Execute the R script with the corresponding files
    Rscript ModDetect.R \
        -f "${BASE_PATH}/${f_files[i]}" \
        -g "${BASE_PATH}/${g_files[i]}" \
        -r "${BASE_PATH}/rep4Uendo.fas" \
        -m 1 \
        -o "${BASE_PATH}/output_ModDetect/${output_files[i]}" \
        -b "${BASE_PATH}/rep4Uendo.bed" \
        -x 1 \
        -y 1 \
        -a
    
    # Add sample information to the output files
    awk -v sample="${samples[i]}" 'BEGIN {FS=OFS=","} {print sample, $0}' "${BASE_PATH}/output_ModDetect/${output_files[i]}" >> "$concat_output"
    
    echo "Completed set $((i + 1))"
done

echo "Concatenation completed. Output saved to: $concat_output"
