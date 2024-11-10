#!/usr/bin/bash
#SBATCH --time=24:00:00
#SBATCH -p normal
#SBATCH -c 32
#SBATCH --mem=32GB

bash bidamp_analyze.sh 20241105_pool1test/recalled/rebasecalling/pass 20241105_pool1test/pool1analysis /oak/stanford/groups/nicolemm/rodell/BIDamplicon/minimap2/20241105/20241105barcode.csv /oak/stanford/groups/nicolemm/rodell/minimap2/GRCh38.fa /oak/stanford/groups/nicolemm/rodell/BIDamplicon/minimap2/20241105/20241105filecombos.csv /oak/stanford/groups/nicolemm/rodell/BIDamplicon/partialpool1.bed

