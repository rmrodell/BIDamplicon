#!/bin/bash

# Define the base directory containing your directories
BASE_DIR="/oak/stanford/groups/nicolemm/rodell/BIDamplicon/20241003/20241003_1058_MN39827_AVY431_84d24661/fastq_pass"
# Define the destination directory for the concatenated files
DEST_DIR="/oak/stanford/groups/nicolemm/rodell/BIDamplicon/minimap2/20241003/"

# Associative array to map barcodes to specific output filenames
declare -A barcode_filenames
barcode_filenames=(
    ["barcode01"]="KD_P101_BS_0902.fq"
    ["barcode02"]="WT_P102_BS_0902.fq"
    ["barcode03"]="WT_P4_BS_0902.fq"
    ["barcode04"]="OE_P3_BS_0902.fq"
    ["barcode05"]="4nMC17_WT_P102_BS_0902.fq"
    ["barcode06"]="40nMC17_WT_P102_BS_0902.fq"
    ["barcode07"]="KD_P101_in_0902.fq"
    ["barcode08"]="WT_P102_in_0902.fq"
    ["barcode09"]="WT_P4_in_0902.fq"
    ["barcode10"]="OE_P3_in_0902.fq"
    ["barcode11"]="4nMC17_WT_P102_in_0902.fq"
    ["barcode12"]="40nMC17_WT_P102_in_0902.fq"

    ["barcode13"]="KD_P101_BS_0906.fq"
    ["barcode14"]="WT_P102_BS_0906.fq"
    ["barcode15"]="WT_P4_BS_0906.fq"
    ["barcode16"]="OE_P3_BS_0906.fq"
    ["barcode17"]="50nMC17_WT_P102_BS_0906.fq"
    ["barcode18"]="500nMC17_WT_P102_BS_0906.fq"
    ["barcode19"]="5uMC17_WT_P102_BS_0906.fq"
    ["barcode20"]="50uMC1&_WT_P102_BS_0906.fq"
    ["barcode21"]="4nMC17_WT_P4_BS_0902.fq"
    ["barcode22"]="40nMC17_WT_P4_BS_0902.fq"
    ["barcode23"]="KD_P101_in_0906.fq"
    ["barcode24"]="WT_P102_in_0906.fq"
    ["barcode25"]="WT_P4_in_0906.fq"
    ["barcode26"]="OE_P3_in_0906.fq"
    ["barcode27"]="50nMC17_WT_P102_in_0906.fq"
    ["barcode28"]="500nMC17_WT_P102_in_0906.fq"
    ["barcode29"]="5uMC17_WT_P102_in_0906.fq"
    ["barcode30"]="50uMC1&_WT_P102_in_0906.fq"
    ["barcode31"]="4nMC17_WT_P4_in_0902.fq"
    ["barcode32"]="40nMC17_WT_P4_in_0902.fq"

    ["barcode33"]="WT1_pLKO_BS_1214.fq"
    ["barcode34"]="KD1_shPUS7_BS_1214.fq"
    ["barcode35"]="WT1_pLKO_in_1214.fq"
    ["barcode36"]="KD1_shPUS7_in_1214.fq"
    ["barcode37"]="WT2_pLKO_BS_1214.fq"
    ["barcode38"]="KD2_shPUS7_BS_1214.fq"
    ["barcode39"]="WT2_pLKO_in_1214.fq"
    ["barcode40"]="KD2_shPUS7_in_1214.fq"
    ["barcode41"]="WT3_pLKO_BS_1214.fq"
    ["barcode42"]="KD3_shPUS7_BS_1214.fq"
    ["barcode43"]="WT3_pLKO_in_1214.fq"
    ["barcode44"]="KD3_shPUS7_in_1214.fq"
)

# List of barcodes
barcodes=("barcode01" "barcode02" "barcode03" "barcode04" "barcode05" "barcode06" "barcode07" "barcode08" 
    "barcode09" "barcode10" "barcode11" "barcode12" "barcode13" "barcode14" "barcode15" "barcode16"
    "barcode17" "barcode18" "barcode19" "barcode20" "barcode21" "barcode22" "barcode23" "barcode24"
    "barcode25" "barcode26" "barcode27" "barcode28" "barcode29" "barcode30" "barcode31" "barcode32"
    "barcode33" "barcode34" "barcode35" "barcode36" "barcode37" "barcode38" "barcode39" "barcode40"
    "barcode41" "barcode42" "barcode43" "barcode44")

# Loop through each barcode directory
for barcode in "${barcodes[@]}"; do
    # Navigate to the directory corresponding to the barcode
    DIR="$BASE_DIR/$barcode"
    
    # Check if the directory exists
    if [ -d "$DIR" ]; then
        echo "Processing $barcode..."

        # Navigate to the directory
        cd "$DIR"

        # Unzip all .gz files
        gunzip *.gz

        # Determine the output filename for this barcode
        output_fastq="${barcode_filenames[$barcode]}"

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
done
