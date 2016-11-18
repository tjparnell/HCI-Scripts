#!/usr/bin/perl

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
foreach my $n (qw(GeneID Length)) {
	# copy the geneID and length columns
	my $i = $Data->find_column($n);
	die "can't find $n column!\n" unless defined $i;
	my $column = $Data->column_values($i);
	$Out->add_column($column);
}

## Identify the count columns
my @indices = ask_user_for_index($Data, 'Enter the count columns  ');


## Calculate the FPKMs
foreach my $index (@indices) {
	
	$Out->add_column( $Data->name($index) . '_FPKM');
	
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

	# calculate and store the fpkms
	$Data->iterate( sub {
		my $row = shift;
		my $length = $row->value(1) || 1;  # how would a zero sneak in here? shouldn't happen 
		my $fpkm = ($row->value($index) * 1_000_000_000) / ($length * $total);
		# put the calculated value in the Output table
		$Out->value($row->row_index, $index, $fpkm);
	} );
}


## Write file
my $s = $Out->save($outfile);
print "Saved file $s\n" if $s;



