#!/bin/bash
#BATCH -J test          #job name
#SBATCH -p fat
#SBATCH -t 72:00:00
#SBATCH -N 1
#SBATCH --cpus-per-task=8
#SBATCH --mem=100G

# ====================== 请在这里修改你的文件名 ======================
R1_RAW=$1 #"input_1.fq.gz"   # 你的 R1 原始文件名
R2_RAW=$2 #"input_2.fq.gz"   # 你的 R2 原始文件名
N_READS=$3             # 下采样条数
SEED=123                 # 随机种子（保证双端配对一致）
outfn=$4 #输出文件路径
# ====================================================================
# 检查参数是否完整
if [ $# -ne 4 ]; then
    echo "使用方法：sbatch $0 R1.fq.gz R2.fq.gz 下采样数 输出目录"
    exit 1
fi
mkdir -p "$outfn"
# 自动提取前缀（自动去掉 _1.fq.gz / _2.fq.gz）
PREFIX=$(basename "$R1_RAW" | sed -e 's/_1.fq.gz$//' -e 's/_1.fastq.gz$//')

# 输出文件名
R1_OUT="${outfn}/${PREFIX}_${N_READS}_1.fq.gz"
R2_OUT="${outfn}/${PREFIX}_${N_READS}_2.fq.gz"

echo "========================================"
echo " PE150 fastq.gz downsampled $N_READS reads for pipeline test"
echo "========================================"
echo "Input R1: $R1_RAW"
echo "Input R2: $R2_RAW"
echo "Output R1: $R1_OUT"
echo "Output R2: $R2_OUT"
echo "Random seed: $SEED"
echo "========================================"

# 开始下采样（同种子，保证配对）
echo "Downsampling R1 ..."
seqkit sample -s $SEED -n $N_READS $R1_RAW -o $R1_OUT

echo "Downsampling R2 ..."
seqkit sample -s $SEED -n $N_READS $R2_RAW -o $R2_OUT

echo -e "\n✅  Downsample completed！"
echo "Outfiles："
echo "  $R1_OUT"
echo "  $R2_OUT"
