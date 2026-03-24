#!/usr/bin/perl  
use strict;  
use warnings;  
 
if (@ARGV != 6) {  
   die "Usage: perl merge_script.pl <file1> <file2> <bc1_col> <bc2_col> <bc3_col> <output_file>\n";  
}  
  
my ($file1_name, $file2_name, $col1, $col2, $col3, $output_name) = @ARGV;

$col1--;
$col2--;
$col3--;
 
 
open my $file1, '<', $file1_name or die "Could not open $file1_name: $!";  
open my $file2, '<', $file2_name or die "Could not open $file2_name: $!";  
open my $output, '>', $output_name or die "Could not open $output_name: $!"; 
my %data;  
 
while (my $line = <$file1>) {  
    chomp $line;
    my @fields = split /\t/, $line;
    if (defined $fields[0]) {   
        $data{$fields[0]} = $fields[1];  
       }
}  
 
while (my $line = <$file2>) {  
      chomp $line;  
      my @fields = split /\t/, $line;   
      print $output "$line\t$data{$fields[$col1]}\t$data{$fields[$col2]}\t$data{$fields[$col3]}\t$data{$fields[$col1]}$data{$fields[$col2]}$data{$fields[$col3]}\n"; 
}
 
close $file1;  
close $file2;  
close $output; 