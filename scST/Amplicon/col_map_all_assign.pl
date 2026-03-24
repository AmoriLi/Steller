#!/usr/bin/perl  
use strict;  
use warnings;  

# 
if (@ARGV != 7) {  
    die "Usage: perl merge_script.pl <file1> <file2> <output_file> <file1_key_column> <file1_new_column> <file2_key_column> <file2_old_column>\n";  
    }  

# 
my ($file1_name, $file2_name, $output_name, $file1_key_col,$file1_new_col, $file2_key_col, $file2_old_col) = @ARGV; 

## 
open my $file1, '<', $file1_name or die "Could not open $file1_name: $!";  
open my $file2, '<', $file2_name or die "Could not open $file2_name: $!";  
open my $output, '>', $output_name or die "Could not open $output_name: $!";
  
$file1_key_col--;  #
$file1_new_col--;
$file2_key_col--;    
$file2_old_col--;

my %data;   

##  
while (my $line = <$file1>) {  
    chomp $line;
    my @fields = split /\t/, $line;
    if (defined $fields[$file1_key_col]) {   
        my $key = $fields[$file1_key_col];
	my $value = $fields[$file1_new_col];
        $data{$key} = $value;   
     }
}  

  
while (my $line = <$file2>) {  
    chomp $line;  
    my @fields = split /\t/, $line;  
    if (defined $fields[$file2_key_col]) {    
        my $id = $fields[$file2_key_col];   
 
        if (exists $data{$id}) {  
            my $new_value = $data{$id};   
            print $output "$line\t$new_value\n";  
        }  else {
           print $output "$line\t$fields[$file2_old_col]\n";
        }
     }  
}

# 
close $file1;  
close $file2;  
close $output;