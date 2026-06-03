#!/bin/bash
#SBATCH -o Top_%j_%a.log
#SBATCH -J Read_check
#SBATCH -p fat,cp2             #use partition intel-sc3
#SBATCH -c 4
#SBATCH --mem=50G

wd=$1
indir=$wd/split_fastq

i=$((SLURM_ARRAY_TASK_ID + 1))
fl=$(sed -n "${i}p" ${wd}/tmp.txt)
echo "Processing sample: $fl"

outdir=${indir}/$fl 
in1=$(ls ${outdir}/*1.fq.gz)

#for read1-polyT beads
#zcat $in1 | awk '{if(FNR%4==2)print $0}' | sort -r | uniq -c > $outdir/r1.uniq.txt

#for formal beads
#Extract sequences after polyT(VN)
zcat $in1 | awk '{if(NR%4==2)print $0}' | \
grep "TTTTTTTTTT" | \
sed 's/.*TTTTTTTT[A,C,G][A,T,C,G]//g' | \
sort | uniq -c | \
sort -r -n -k 1 | \
head -n 50 > $outdir/r1.uniq.sort.top50_after_polyTVN.txt

echo "Finished: $outdir/r1.uniq.sort.top50_after_polyTVN.txt"

