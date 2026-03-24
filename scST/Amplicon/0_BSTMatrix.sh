#!/bin/bash
#SBATCH -J BSTMatrix            #job name
#SBATCH -p cp2
#SBATCH -t 72:00:00
#SBATCH -N 1
#SBATCH --cpus-per-task=8
#SBATCH --mem=150G



wkdir=$1 
refdir=~/reference/star/mouse/P045_genocode_vM23/P045_mcherry_BC
gtf=$refdir/custom_gencode.vM23.annotation.gtf
feature=$refdir/features.tsv

cd $wkdir
outdir=$wkdir/BST
#indir=$wkdir/$fl/amplicon/fastq
fqdir=$wkdir/fastq  
fdir=`dirname $wkdir`
#outdir=$wkdir/BST
if [ ! -d $outdir ]; then
mkdir -p $outdir
fi

fq1=`ls ${fqdir}/*1.fq.gz`
fq2=`ls ${fqdir}/*2.fq.gz`
flu=`ls ${fdir}/*FL.tif`
#he=`ls ${fdir}/*HE.tif`
mask=`ls ${fdir}/*-final_signal.txt`

#cp ~/NPC_project/BMK/cs_design/240110_batch1_test2/config.txt $outdir/config.txt
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
#echo "HE he" >> $outdir/config.txt
#echo "CellSplit True" >> $outdir/config.txt
echo "fluorescence $flu" >> $outdir/config.txt
echo "fluorescence_channl 2" >> $outdir/config.txt



BSTMatrix -c $outdir/config.txt -s 1,2,3,4,5