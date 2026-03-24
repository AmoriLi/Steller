#!/usr/bin/perl 
use strict;  
use warnings;  

if (@ARGV != 3) {  
    die "Usage: perl seq_stat.pl <input_file> <output_file> <total_reads>\n";  
    }  

$| = 1;

my $infile = shift @ARGV;  
my $outfile = shift @ARGV;
my $total_reads =  shift @ARGV;

&get_umi_gene_and_stat($infile, $outfile, $total_reads);

sub get_umi_gene_and_stat()
{
    my ($infile, $outfile, $total_reads) = @_ ;

    open (IN, $infile) || die "$infile, $!\n" ;
    my $uniq = 0 ;
    my $total = 0 ;
    my @simu_arr = ();
    my $index = 0 ;
    while(<IN>){
        chomp ;
        next if (m/^\#/);
        $index++ ;
        my ($bc, $umi, $gene, $num) = split ;
        $uniq++ ;
        $total += $num ;
        push @simu_arr, ($index)x$num ;
    }
    close(IN);
    &shuffle(\@simu_arr);

    my $percent = 0 ;
    my @uniqs = ();
    my $simu_uniq = 0 ;
    my %hstat = ();
    for (my $i=0; $i<@simu_arr; $i++){
        if (!defined $uniqs[$simu_arr[$i]]){
            $uniqs[$simu_arr[$i]] = 1 ;
            $simu_uniq++ ;
        }
        if ($i+1 >= $total*$percent/100){
            my $uniq = $simu_uniq ;
            my $simu_total = $i + 1 ;
            $hstat{$percent} = [$uniq, $simu_total, int((1-$uniq/$simu_total)*10000+0.5)/100] ;
            $percent += 5 ;
        }
    }

    open (OUT, ">$outfile") || die "$outfile, $!\n" ;
    print OUT "#Percent(%)\tUniq\tTotal\tSeqSaturateion(%)\n" ;
    $total_reads ||= $total ;
    for my $perc (sort {$a<=>$b} keys %hstat){
        my ($simu_uniq, $simu_total, $simu_saturat) = @{$hstat{$perc}} ;
        my $tr_uniq = int($simu_uniq/$total*$total_reads) ;
        my $tr_total = int($simu_total/$total*$total_reads) ;
        print OUT "$perc\t$tr_uniq\t$tr_total\t$simu_saturat\n" ;
    }
    close(OUT);

    return ;
}

sub shuffle()
{
    my ($asimu) = @_ ;

    my $num = scalar(@$asimu);
    for (my $i=0; $i<@$asimu; $i++){
        my $index = int(rand($num));
        ($asimu->[$i], $asimu->[$index]) = ($asimu->[$index], $asimu->[$i]);
    }

    return ;
}