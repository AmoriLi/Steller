#!/bin/bash
#SBATCH -J fastqc
#SBATCH -p cp2
#SBATCH -t 72:00:00
#SBATCH -N 1
#SBATCH --cpus-per-task=8
#SBATCH --mem=100G

script="~/NPC_project/P045_whitelist/pipeline"
wd=$1 #e.g. ~/NPC_project/P045_whitelist/w2
cd $wd
fqdir=$wd/fastq
outdir=$wd/fastp
if [ ! -d "$outdir" ]; then
	mkdir -p "$outdir"
fi

#i=`expr ${SLURM_ARRAY_TASK_ID} + 1`
for file in `ls $fqdir/*.fastq.gz | sed 's/_R[1,2]_001.fastq.gz//g'| uniq`
do 
sample=$file
echo $sample
fq1=${sample}_R1_001.fastq.gz
fq2=${sample}_R2_001.fastq.gz


#merge PE150
#seqtk mergepe "$fq1" "$fq2" | gzip > merged.fastq.gz
fastp -i "$fq1" -I "$fq2" --merge --merged_out $outdir/merged.fastq --out1 $outdir/unmerged_r1.fastq --out2 $outdir/unmerged_r2.fastq

done
