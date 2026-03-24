#!/bin/bash
#BATCH -J starcode         #job name
#SBATCH -p fat,cp2
#SBATCH -t 72:00:00
#SBATCH -N 1
#SBATCH --cpus-per-task=4
#SBATCH --mem=100G
wd=$1
starcode -s -d 5 -t 4 --print-clusters $wd/read_filtered.clonebc.txt > $wd/read_filtered.clonebc.collapsed.txt