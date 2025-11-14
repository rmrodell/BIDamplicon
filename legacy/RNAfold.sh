#!/bin/bash

ml biology bedtools viennarna
ml R

DIR=$(pwd)

for bed in "$DIR"/*.bed; do
    
    sample=$(basename "${bed}" .bed)  # Get the base name without the extension

    echo "Processing $sample"

    mkdir -p $sample # make directory to sort sample into

    bedtools getfasta -fi /home/users/rodell/pool1/pool1_prettyplease.fasta -bed $bed -fo $sample/$sample.fa # create fasta

    (cd $sample # change directory for processing
    
    RNAfold -p -o --MEA $sample.fa # fold everything in the fasta

    /scratch/users/rodell/RNAfold/RNAFold.R # pull out the important values
    )
done