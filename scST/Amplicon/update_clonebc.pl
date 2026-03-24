#!/usr/bin/perl


use strict;
use warnings;

if (@ARGV != 3) {
    die "Usage: perl script.pl <text1> <text2> <outfn>\n";  
}

my ($text1, $text2, $outfn) = @ARGV; 

open my $infile1, '<', $text1 or die "Could not open $text1: $!";  

my %mapping;

while (my $line = <$infile1>) {
    chomp $line;
    my @fields = split(/\t/, $line);  

    my $col1 = $fields[0];
    my @listA = split(/,/, $fields[2]);
    foreach my $item (@listA) {
        $mapping{$item} = $col1;
    }
}

close($infile1);

open my $infile2, '<', $text2 or die "Could not open $text2: $!"; 

open my $outfile, '>', $outfn or die "Could not open collasped_$text2: $!";

while (my $line = <$infile2>) {
    chomp $line;
    my @fields = split(/\t/, $line);  
    if (exists $mapping{$fields[1]}) {
        $fields[1] = $mapping{$fields[1]};
    }
    print $outfile join("\t", @fields), "\n";
    
}

close($infile2);
close($outfile);

print "Processing complete. Results are in $outfile.\n";