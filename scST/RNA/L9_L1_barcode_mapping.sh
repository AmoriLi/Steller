#!/bin/bash
#BATCH -J barcode_map     #job name
#SBATCH -p fat
#SBATCH -t 72:00:00
#SBATCH -N 1
#SBATCH --cpus-per-task=16
#SBATCH --mem=50G
#SBATCH -a 0-1
wd=~/NPC_project/BMK/tissue/brain/P5
script=~/NPC_project/BMK/pipeline/high_level_level1_mapping.pl
i=`expr ${SLURM_ARRAY_TASK_ID} + 1`
fl=$(echo `sed -n "${i}p" tmp.txt`)
echo $fl
infile1=${wd}/${fl}/RNA/BST/05.AllheStat/level_matrix/level_1/barcodes_cluster.tsv.gz
infile2=${wd}/${fl}/RNA/BST/05.AllheStat/level_matrix/level_9/barcodes_cluster.tsv.gz
outfile=~/NPC_project/BMK/tissue/brain/P5/Merged_2512/2_region/L9/${fl}_L9_map2_L1.txt

perl $script $infile1 $infile2 L9 > $outfile 