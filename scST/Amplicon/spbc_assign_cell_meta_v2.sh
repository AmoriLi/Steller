#!/bin/bash 

##assign L1, cell label and position
#assign L1_label and position
wd=`dirname $1` #~/NPC_project/BMK/tissue/brain/E15.5/2024.5.7/S3000-A4
tmpdir=$2 #$wd/amplicon/result/CloneBC_v4
l1_pre=$wd/RNA/BST/05.AllheStat/level_matrix/level_1
cell_pre=$wd/RNA/BST/07.CellSplit
RNA_spBC=$wd/RNA/BST/RNA/barcode.tsv
script=~/NPC_project/BMK/pipeline

awk '{print NR "\t" $0}' $RNA_spBC > $tmpdir/barcode_id.tsv
zcat $l1_pre/barcodes_cluster.tsv.gz > $tmpdir/barcodes_cluster.tsv
awk -F'\t' -v OFS='\t' 'NR==FNR{a[$1]=$2;next} {if($1 in a) {print $1,a[$1],"L1_" $2}}' $tmpdir/barcode_id.tsv $tmpdir/barcodes_cluster.tsv > $tmpdir/barcodes_id_level1.tsv
zcat $l1_pre/barcodes_pos.tsv.gz > $tmpdir/barcodes_pos.tsv
awk -F'\t' -v OFS='\t' 'NR==FNR{a[$1]=$2"\t"$3;next} {if($3 in a) {print $0,a[$3]}}' $tmpdir/barcodes_pos.tsv $tmpdir/barcodes_id_level1.tsv  > $tmpdir/barcodes_id_level1pos.tsv
bash $script/inner_merge_file.sh $cell_pre/cell_split_result/all_barcode_num.txt $tmpdir/barcodes_id_level1pos.tsv 1 2 3 > $tmpdir/barcodes_id_level1pos_cell.tsv
awk -F'\t' -v OFS='\t' 'NR==FNR{a[$1]=$2"\t"$3;next} {if($6 in a) {print $0,a[$6]}}' $cell_pre/mtx/cells_center.txt $tmpdir/barcodes_id_level1pos_cell.tsv  > $tmpdir/barcodes_id_level1pos_cellpos.tsv
awk -F'\t' -v OFS='\t' 'NR==FNR{a[$2]=$6"\t"$7"\t"$8;next} {if($2 in a) {print $0,a[$2]}}' $tmpdir/barcodes_id_level1pos_cellpos.tsv $tmpdir/spbc2.umi2.clonebc > $tmpdir/spbc2.umi2.clonebc.cellmeta.tsv