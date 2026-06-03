#!/bin/bash
#SBATCH -o fastp_%j.log         
#SBATCH -p fat,cp2             #use partition intel-sc3
#SBATCH -c 8
#SBATCH --mem=50G
wd=$1
inpath=${wd}/fastq
outdir=${wd}/fastp

if [ ! -d $outdir ];then
mkdir -p $outdir
fi

fq1=`ls ${inpath}/*_1.fq.gz`

sample=$(basename $fq1 _1.fq.gz)
out1=${outdir}/${sample}_Q25_clean_R1.fq.gz

fastp -i $fq1 \
-o $out1 \
-q 25 \
-w 8 \
-h ${outdir}/${sample}_clean_report.html \
-j ${outdir}/${sample}_clean_report.json #\

echo "Clean fastq generated: $out1"
echo "Fastp finished successfully!"
