#!/usr/bin/perl

# Timothy J. Parnell, PhD
# Huntsman Cancer Institute
# University of Utah
# Salt Lake City, UT 84112
#  
# This package is free software; you can redistribute it and/or modify
# it under the terms of the Artistic License 2.0.  
# 
# Updated versions of this file may be found in the repository
# https://github.com/tjparnell/HCI-Scripts/

use strict;
use Bio::ToolBox::Data;
use Bio::ToolBox::utility;

unless (@ARGV) {
	print <<HELP;
A script to calculate fpkm values from a featureCounts table.

The program expects a "geneid" and "length" column and will keep 
those columns. It will interactively ask the user for the counts 
columns. FPKM values are written to a second file.

Usage: $0 <counts.txt> <fpkm.txt>

HELP
	exit;
}


## Open input file
my $file = shift @ARGV;
my $Data = Bio::ToolBox::Data->new(file => $file) or 
	die "unable to open input file '$file'\n";
printf " Opened '$file' with %d columns and %d data rows\n", 
	$Data->number_columns, $Data->last_row;

## Prepare output file
my $outfile = shift @ARGV || $Data->path . $Data->basename . '_FPKM.txt';
my $Out = Bio::ToolBox::Data->new();

## Find columns
my $gene_i = $Data->find_column('GeneID') || $Data->name_column;
die " can't find GeneID or Name column\n" unless defined $gene_i;
my $length_i = $Data->find_column('Length');
die " can't find Length column\n" unless defined $length_i;

## Copy columns
foreach my $i ($gene_i, $length_i) {
	# copy the geneID and length columns
	my $column = $Data->column_values($i);
	$Out->add_column($column);
}

## Identify the count columns
my @indices = ask_user_for_index($Data, 'Enter the count columns  ');


## Calculate the FPKMs
foreach my $index (@indices) {
	
	my $new_index = $Out->add_column( $Data->name($index) . '_FPKM');
	
	# determine total number of counts for this dataset
	# we're going line by line to check for potential non-integers
	my $total = 0;
	$Data->iterate( sub {
		my $row = shift;
		my $v = $row->value($index);
		$total += $v if $v =~ /^\d+$/;
	} );
	die sprintf("The sum of column %s is zero! Check your file!\n", 
		$Data->name($index)) unless $total;
	printf "  processing %s, total count $total\n", $Data->name($index);
	
	# calculate and store the fpkms
	$Data->iterate( sub {
		my $row = shift;
		my $length = $row->value($length_i) || 1;  # how would a zero sneak in here? shouldn't happen 
		my $fpkm = ($row->value($index) * 1_000_000_000) / ($length * $total);
		# put the calculated value in the Output table
		$Out->value($row->row_index, $new_index, $fpkm);
	} );
}


## Write file
my $s = $Out->save($outfile);
print "Saved file $s\n" if $s;



