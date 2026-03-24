#!/bin/bash
#BATCH -J SeekSoul        #job name
#SBATCH -p fat,cp2
#SBATCH -t 72:00:00
#SBATCH -N 1
#SBATCH --cpus-per-task=8
#SBATCH --mem=100G

fqdir=$1
refdir=~/reference/star/mouse/P045_genocode_vM23/P045_mcherry_BC
gtf=$refdir/custom_gencode.vM23.annotation.gtf
sn=$2
outdir=$1/result

fq1=`ls ${fqdir}/*expression/*R1_001.fastq.gz`
fq2=`ls ${fqdir}/*expression/*R2_001.fastq.gz`
spfq1=`ls ${fqdir}/*spatial/*R1_001.fastq.gz`
spfq2=`ls ${fqdir}/*spatial/*R2_001.fastq.gz`
hdmifq=`ls ${fqdir}/*expression/*.fq.gz`
chipID=$(basename "$hdmifq" .fq.gz)
dapi=`ls ${fqdir}/*expression/${chipID}.tif`
he=`ls ${fqdir}/*expression/*_HE.tif`

echo $fq1
echo $fq2
echo $spfq1
echo $spfq2
echo $sn
echo $chipID

seekspacetools run \
--fq1 $fq1 \
--fq2 $fq2 \
--spatialfq1 $spfq1 \
--spatialfq2 $spfq2 \
--hdmifq $hdmifq \
--samplename $sn \
--outdir $outdir \
--genomeDir $refdir \
--gtf $gtf \
--chemistry DDVS \
--core 8 \
--include-introns \
--min_umi 1 \
--chip_id $chipID \
--DAPI $dapi \
--HE $he