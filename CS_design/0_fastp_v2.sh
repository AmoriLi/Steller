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

fq1=$(ls ${inpath}/*_1.fq.gz)
fq2=$(ls ${inpath}/*_2.fq.gz)

sample=$(basename $fq1 _1.fq.gz)

fastp -i $fq1 -I $fq2 \
-l 140 \
-Y 30 \
-o $outdir/${sample}_clean.r1.fq.gz \
-O $outdir/${sample}_clean.r2.fq.gz \
-q 25 \
-u 30 \
-w 8 \
-m \
--merged_out ${outdir}/${sample}_merge_clean.fastq.gz \
-h ${outdir}/${sample}_clean.html \
-j ${outdir}/${sample}_clean.json
 
#-m --merged_out ${outdir}/${sample}_merge_clean.fastq.gz \
#--failed_out ${outdir}/${sample}_failed.fastq.gz \
#--include_unmerged \
#--overlap_len_require 6 \
#--overlap_diff_percent_limit 20 \
#--detect_adapter_for_pe #\
#-5 \
#-r \
#-l 50 \
#-n 5\
#-y 

