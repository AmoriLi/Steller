#!/bin/bash
#SBATCH -J BSTMatrix            #job name
#SBATCH -p fat
#SBATCH -t 72:00:00
#SBATCH -N 1
#SBATCH --cpus-per-task=16
#SBATCH --mem=200G

wkdir=$1 #~/NPC_project/BMK/tissue/brain/P5/S3000/RNA
refdir=$2 #~/reference/star/mouse/P045_genocode_vM23/P045_mcherry_BC
gtf=$refdir/custom_gencode.vM23.annotation.gtf
feature=$refdir/features.tsv
indir=$wkdir/fastq
fdir=`dirname $wkdir`
outdir=$wkdir/BST
if [ ! -d $outdir ]; then
mkdir -p $outdir
fi

fq1=`ls ${indir}/*1.fq.gz`
fq2=`ls ${indir}/*2.fq.gz`
flu=`ls ${fdir}/*FL.tif`
#he=`ls ${fdir}/*HE.tif`
mask=`ls ${fdir}/*-final_signal.txt`

echo "BCType V2" >> $outdir/config.txt
echo "BCThreads 8" >> $outdir/config.txt
echo "Sjdboverhang 100" >> $outdir/config.txt 
echo "STARThreads 8" >> $outdir/config.txt
echo "FQ1 $fq1" >> $outdir/config.txt
echo "FQ2 $fq2" >> $outdir/config.txt
echo "INDEX $refdir/" >> $outdir/config.txt
echo "GFF $refdir/custom_gencode.vM23.annotation.gtf" >> $outdir/config.txt
echo "FEATURE $feature" >> $outdir/config.txt
echo "OUTDIR $outdir" >> $outdir/config.txt
echo "PREFIX RNA" >> $outdir/config.txt
echo "FLU $mask" >> $outdir/config.txt
#echo "HE $he" >> $outdir/config.txt
echo "CellSplit True" >> $outdir/config.txt
echo "fluorescence $flu" >> $outdir/config.txt
echo "fluorescence_channl 2" >> $outdir/config.txt
#specify cell split parameters
echo "diameter [28]" >> $outdir/config.txt
echo "cell_edge 1" >> $outdir/config.txt #without expand to cytoplasma
echo "min_detect_value 50" >> $outdir/config.txt
echo "flow_threshold 0.2" >> $outdir/config.txt


BSTMatrix -c $outdir/config.txt -s 1,2,3,4,5,7
#BSTMatrix -c $outdir/config.txt -s 7
