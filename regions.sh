# this version does NOT pay attention to strandedness
BED=$1 # first argument is bed file with sites of interest

ml biology bedtools

echo "Examining $BED"

# intersect with 5'UTR
num_5utr=$(bedtools intersect -a "$BED" -b /home/users/rodell/regions/hg38_5utr.bed -u | wc -l)
echo "Number of sites in 5' UTR: $num_5utr"

# interest with 3' UTR
num_3utr=$(bedtools intersect -a "$BED" -b /home/users/rodell/regions/hg38_3utr.bed -u | wc -l)
echo "Number of sites in 3' UTR: $num_3utr"

# interest with exons
num_exons=$(bedtools intersect -a "$BED" -b /home/users/rodell/regions/hg38_codingexons.bed -u | wc -l)
echo "Number of sites in exons: $num_exons"

# interest with 3' UTR
num_introns=$(bedtools intersect -a "$BED" -b /home/users/rodell/regions/hg38_introns.bed -u | wc -l)
echo "Number of sites in introns: $num_introns"