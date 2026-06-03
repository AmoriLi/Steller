#!/bin/bash
#SBATCH -o STAR_%j_%a.log
#SBATCH -J STAR            #job name
#SBATCH -p fat,cp2             #use partition intel-sc3
#SBATCH -c 10
#SBATCH --mem=50G
wd=$1
indir=$wd/split_fastq
i=$((SLURM_ARRAY_TASK_ID + 1))
fl=$(sed -n "${i}p" ${wd}/tmp.txt)
echo "Processing sample: $fl"

r1=$indir/$fl/${fl}.r1.fq.gz
r2=$indir/$fl/${fl}.r2.fq.gz
outdir=$indir/$fl/STAR

if [ ! -d $outdir ]; then
mkdir -p $outdir
fi

#reference genome
refdir=$2 #star/mouse/P045_genocode_vM23/P045_mcherry_star_index
gtf=$refdir/custom_gencode.vM23.annotation.gtf

# ============ STAR mapping ==============
STAR --runMode alignReads \
--outSAMtype BAM SortedByCoordinate \
--runThreadN 10 \
--genomeDir $refdir \
--readFilesIn $r1 $r2 \
--outFileNamePrefix $outdir/ \
--readFilesCommand zcat \
--outFilterMatchNminOverLread 0.33 \
--outFilterScoreMinOverLread 0.33 
#--sjdbOverhang 100 \
#--outSJfilterReads Unique

# ============== Feature counts =======================
featureCounts -p \
    -a $gtf \
    -T 10 \
    -t gene \
    -g gene_id \
    -o $outdir/count.txt \
    $outdir/Aligned.sortedByCoord.out.bam 

echo "Done: $fl"
