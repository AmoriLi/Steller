#!/bin/bash
#BATCH -J bam_process          #job name
#SBATCH -p cp2
#SBATCH -t 72:00:00
#SBATCH -N 1
#SBATCH --cpus-per-task=8
#SBATCH --mem=100G

#
start_time=$(date "+%s")

pt='[ATCG]{4}CT[ATCG]{4}AC[ATCG]{4}TC[ATCG]{4}GT[ATCG]{4}TG[ATCG]{4}CA[ATCG]{4}'
wd=$1 #"~/NPC_project/BMK/tissue/brain/E15.5/2024.5.7/S3000-A4/amplicon"
test=$2 #T or F
#collapse=$3 #if do starcode
#collapse_type=$4 # select which group to do starcode collapse
indir=$wd/fastq
indir1=$wd/BST/01.fastq2BcUmi
indir2=$wd/BST/02.Umi2Gene
bcref="~/NPC_project/BMK/pipeline/BSTMatrix_v2.4.e/ref_file/V2bcs.fa"
whitelist=~/NPC_project/P045_whitelist/P045/whitelist.tsv
script="~/NPC_project/BMK/pipeline/amplicon"
#outdir=$wd/result/CloneBC/1_preprocess
tmpdir=$wd/tmp
if [ ! -d $tmpdir ]; then
    mkdir -p $tmpdir
fi

if [ ! -d $outdir ]; then
    mkdir -p $outdir
fi
n=50000
fq2=`ls ${indir}/*2.fq.gz`

if [ "$test" = "T" ]; then
   seqtk sample $fq2 $n | gzip -c > $tmpdir/test100w.r2.fq.gz
   ##add id to fa
   zcat $tmpdir/test100w.r2.fq.gz | awk 'NR%4==2{print int(NR/4)"\t"$0}' > $tmpdir/read.id.r2.fq
   ##filter mcherry bam id
   samtools view $indir2/RNAAligned.out.bam | awk '{if ($3=="mcherry") print $1}' | uniq > $tmpdir/mcherry.bam.id
else
   #zcat $fq2 | awk 'NR%4==2{print int(NR/4)"\t"$0}' > $tmpdir/read.id.r2.fq
   samtools view $indir2/RNAAligned.out.bam | awk '{if ($3=="mcherry") print $1"\t"$10}' | uniq > $tmpdir/mcherry.read.id.r2.fq
fi
##filter fastq
#store total fastq read number into stat.txt
n=$(zcat "$fq2" | wc -l)
n=$((n / 4))
printf "%s\t%d\n" "fastq" "$n" > $tmpdir/rd.stat.txt 

#store mcherry aligned read number into stat.txt
n=$(cat $tmpdir/mcherry.read.id.r2.fq | wc -l)
printf "%s\t%d\n" "mcherry" "$n" >> $tmpdir/rd.stat.txt
##final readid, spbc, umi2 file
perl $script/col_filter_by_file.pl $tmpdir/mcherry.read.id.r2.fq $indir1/RNA.umi $tmpdir/RNA.umi2 1 1
awk '{print $1"\t"$2"\t"$2"_"$5"\t"$5}' $tmpdir/RNA.umi2 > $tmpdir/readid.spbc.umi.txt 
perl $script/col_map_all_assign.pl $indir1/RNA.umi_cor.info $tmpdir/readid.spbc.umi.txt $tmpdir/readid.spbc.umi2 1 4 3 4

#recommended!: qc1 by umi nread>1
python $script/umi_nread_check_and_filter.py $tmpdir 1
#read number with confident umi2 read
n=$(cat $tmpdir/readid.spbc.umi2.nreadcut1.tsv | wc -l)
printf "%s\t%d\n" "umi_passed_nread" "$n" >> $tmpdir/rd.stat.txt

#map spbc sequence
awk '{split($2,a,"-");print $0"\t"a[1]"\t"a[2]"\t"a[3]}' $tmpdir/readid.spbc.umi2.nreadcut1.tsv > $tmpdir/readid.splitbc.umi2
#generate spbc, sequence file
awk -F '>' '{if (NF>1) {bc=$2;next} else {print bc"\t"$0}}' $bcref > $tmpdir/V2bcs.txt
perl $script/get_spbc_seq.pl $tmpdir/V2bcs.txt $tmpdir/readid.splitbc.umi2 8 9 10 $tmpdir/readid.spbc2.umi2

###extract clonebc
#add readsequence
perl $script/col_map_filter_assign.pl $tmpdir/mcherry.read.id.r2.fq $tmpdir/readid.spbc2.umi2 $tmpdir/readid.read.spbc2.umi2 1 2 1

perl $script/larry_pattern_sub.pl $tmpdir/readid.read.spbc2.umi2 15 $tmpdir/spbc2.umi2.clonebc

n=$(cat $tmpdir/spbc2.umi2.clonebc | wc -l)
printf "%s\t%d\n" "clonebc" "$n" >> $tmpdir/rd.stat.txt

# sequencing saturation based on valid spbc_umi_clonebc nreads
### spbc2.umi2.clonebc can be used to filter with whitelist CloneBC and get final valid spbc_umi_clonebc data, then calculate sequencing saturation curve
perl $script/col_filter_by_file.pl $whitelist $tmpdir/spbc2.umi2.clonebc $tmpdir/spbc2.umi2.clonebc2 1 16 

awk '
{
    key = $2 "\t" $5 "\t" $14
    count[key]++
}
END {
    for (key in count) {
        print key "\t" count[key]
    }
}' $tmpdir/spbc2.umi2.clonebc > $tmpdir/spbc2.umi2.clonebc2.nread

perl $script/seq_sat.pl $tmpdir/spbc2.umi2.clonebc2.nread $tmpdir/sequence_saturation.stat $n

bash $script/spbc_assign_cell_meta_v2.sh $wd $tmpdir
n=$(cat $tmpdir/spbc2.umi2.clonebc.cellmeta.tsv | wc -l)
printf "%s\t%d\n" "cell" "$n" >> $tmpdir/rd.stat.txt

### the parameter can be manually adjust by running this script interactively.
python $script/nread_check_and_filter_v3.py $tmpdir 5 10 #5
n=$(cat $tmpdir/umi2.clonebc.cellmeta.nreadcut.tsv | wc -l)
printf "%s\t%d\n" "cell_nread" "$n" >> $tmpdir/rd.stat.txt

python $script/clonebc_count_generate_v2.py $tmpdir

starcode -s --print-clusters $tmpdir/clonebc_count.txt > $tmpdir/clonebc_count_collapsed.txt

perl $script/update_clonebc.pl $tmpdir/clonebc_count_collapsed.txt $tmpdir/umi2.clonebc.cellmeta.nreadcut.tsv $tmpdir/collapsed_clonebc_umi_cell.txt
# 
end_time=$(date "+%s")

# 
duration=$(( (end_time - start_time) / 60 ))
echo "CloneBC preprocess duration: $duration min"
