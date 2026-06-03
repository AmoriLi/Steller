#!/bin/bash
#SBATCH -o a%j_%a.log
#SBATCH -J split_fastq            #job name
#SBATCH -p fat,cp2             #use partition intel-sc3
#SBATCH -c 8
#SBATCH --mem=50G
wd=$1
indir=$wd/fastq
in1=`ls ${indir}/*1.fq.gz`
in2=`ls ${indir}/*2.fq.gz`
outdir=$wd/split_fastq

#all probe samplenames in split_fastq directory
cd $outdir
ls -d */ | sed 's/\///g' | sort | uniq > $wd/tmp.txt
cd -

i=$((SLURM_ARRAY_TASK_ID + 1)) #`expr ${SLURM_ARRAY_TASK_ID} + 1`
fl=$(sed -n "${i}p" $wd/tmp.txt)
echo "Processing sample: $fl"

sample_dir=${outdir}/$fl
readname_file=${sample_dir}/${fl}_readname.txt
r1=${sample_dir}/${fl}.r1.fq
r2=${sample_dir}/${fl}.r2.fq

if [ ! -d $sample_dir ]; then
mkdir -p $sample_dir
fi

#extract fastq records
seqkit grep -f $readname_file $in1 > $r1
seqkit grep -f $readname_file $in2 > $r2

#compress
gzip -f $r1
gzip -f $r2

echo "Done: $fl"
