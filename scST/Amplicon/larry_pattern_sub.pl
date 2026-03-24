#!/usr/bin/perl

use strict;
use warnings;

if (@ARGV != 3) {  
     die "Usage: perl script.pl <infile> <column> <outfile>\n";  
         }
  
my ($infn, $col, $outfn) = @ARGV; 

my $pattern = '[ATCG]{4}CT[ATCG]{4}AC[ATCG]{4}TC[ATCG]{4}GT[ATCG]{4}TG[ATCG]{4}CA[ATCG]{4}';
 
open my $infile, '<', $infn or die "Could not open $infn: $!";  
open my $outfile, '>', $outfn or die "Could not open $outfn: $!";  
 
$col--;   

while (my $line = <$infile>) {
    chomp $line;
    my @fields = split(/\t/, $line); 
    if (defined $fields[$col]) {
      if ($fields[$col] =~ /$pattern/) {
          print $outfile "$line\t$&\n";                                
       }
    }
}

close($infile);
close($outfile);