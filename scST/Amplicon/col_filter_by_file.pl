#!/usr/bin/perl  
use strict;  
use warnings;  
##  
if (@ARGV != 5) {  
    die "Usage: perl merge_script.pl <file1> <file2> <output_file> <file1_key_column> <file2_key_column>\n";  
        }  
#  
my ($file1_name, $file2_name, $output_name, $file1_col, $file2_col) = @ARGV; 
##
open my $file1, '<', $file1_name or die "Could not open $file1_name: $!";  
open my $file2, '<', $file2_name or die "Could not open $file2_name: $!";  
open my $output, '>', $output_name or die "Could not open $output_name: $!";
# index start from 0  
$file1_col--;
$file2_col--;    
my %data;  # key-value  
## file1  
while (my $line = <$file1>) {  
   chomp $line;
   my @fields = split /\t/, $line;
   if (defined $fields[$file1_col]) {   
      $data{$fields[$file1_col]} = 1;  #   
     }
}  
# file2 only matched file1   
while (my $line = <$file2>) {  
    chomp $line;  
    my @fields = split /\t/, $line;  #  
    if (defined $fields[$file2_col]) {  # 
       if (exists $data{$fields[$file2_col]}) { 
            print $output "$line\n"; 
        }  
    }  
}

#   
close $file1;  
close $file2;  
close $output; 