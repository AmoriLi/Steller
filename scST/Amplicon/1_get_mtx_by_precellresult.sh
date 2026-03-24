#!/bin/bash
#SBATCH -J cellsplit_mtx            #job name
#SBATCH -p fat
#SBATCH -t 72:00:00
#SBATCH -N 1
#SBATCH --cpus-per-task=8
#SBATCH --mem=100G

#2024-12-02
#get amplicon cell/nuceli split matrix from RNA cell/nuclei split npy file

wkdir=$1 #~/NPC_project/BMK/tissue/brain/P5/S3000_50cs3/amplicon
cd $wkdir

fdir=`dirname $wkdir`

wd=$fdir/RNA #~/NPC_project/BMK/tissue/brain/P5/S3000_50cs3/RNA
Bin=~/software/BSTMatrix_v2.4.e
outdir=$wkdir/BST/07.CellSplit
npyfile=$wd/BST/07.CellSplit/cell_split_result/cells.npy

if [ ! -d $outdir ]; then
mkdir -p $outdir/mtx
fi

python $Bin/cell_split/get_mtx.py -i $wkdir/BST/05.AllheStat/level_matrix/level_1  -c $wd/BST/07.CellSplit/cell_split_result/all_barcode_num.txt -o $outdir/mtx

python $Bin/cell_split/find_cells_loc.py --cells_path $npyfile --out_path $outdir/mtx/

#perl $Bin/cell_split/generate_bc_pos.pl $outdir/mtx/barcodes.tsv.gz $outdir/mtx/cells_center.txt $outdir/mtx/barcodes_pos.tsv && gzip -f $outdir/mtx/barcodes_pos.tsv