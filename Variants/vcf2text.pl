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

unless (@ARGV) {
	print <<HELP;
Export VCF to tab delimited file with defined sample attributes.
Usage:
 $0 <file.vcf> <output.txt> <tag> <tag> ...
HELP
	exit;
}

# parameters
my $infile = shift @ARGV;
my $outfile = shift @ARGV or die "no output file provided!\n";
my @tags = @ARGV or die "no sample FORMAT tags provided!\n";

# input file stream
my $In = Bio::ToolBox::Data->new(file => $infile, stream => 1) or 
	die "unable to open file $infile!\n";
die "$infile is not VCF format!\n" unless $In->vcf;

# output file
my $Out = Bio::ToolBox::Data->new(
	stream => 1, 
	columns => [qw(Chrom Pos ID Allele)],
	out => $outfile,
) or die "unable to make output data object!\n";
my $number = $In->number_columns - 1;
for my $i (9 .. $number) {
	foreach my $t (@tags) {
		$Out->add_column( sprintf("%s_%s", $In->name($i), $t) );
	}
}	

# add data
while (my $row = $In->next_row) {
	# get chromo, pos, alt allele
	my @data = ($row->seq_id, $row->start, $row->value(2), $row->value(4)); 
	my $attrib = $row->attributes;
	
	# walk through samples and collect attributes
	for my $i (9 .. $number) {
		foreach my $t (@tags) {
			push @data, $attrib->{$i}{$t} || 0;
		}
	}	
	
	# write
	$Out->write_row(\@data);
}

# finish
$In->close_fh;
$Out->close_fh;


