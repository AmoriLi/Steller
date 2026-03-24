#!/usr/bin/perl
use strict;
use warnings;

if(@ARGV!=3){
	print "perl $0 <BSTViewer_project/level_matrix/level_1/barcodes_cluster.tsv.gz> <BSTViewer_project/level_matrix/level_4/barcodes_cluster.tsv.gz> <L4>\n";
	exit;
}
open(IN,"gunzip -c $ARGV[0]|");
my %level1;
while(<IN>){
	chomp;
	my @a=split(/\t/,$_);
	$level1{$a[0]}="L1_".($a[1]-1);
}

open(IN2,"gunzip -c $ARGV[1]|");
while(<IN2>){
	chomp;
	my @a=split(/\t/,$_);
	my $L1=$level1{$a[0]};
	print "$L1\t$ARGV[2]_".($a[1]-1),"\n";
}
