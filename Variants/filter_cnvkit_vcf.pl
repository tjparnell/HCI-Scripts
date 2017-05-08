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
use Bio::ToolBox::Data::Stream;

unless (scalar(@ARGV) >= 2) {
	print <<HELP;

Filter a VCF file generated from CNVkit of copy number alterations by it's 
log fold change Info attribute FOLD_CHANGE_LOG. Keeps all variants whose 
absolute value is greater than or equal to the provided threshold.

Usage: $0 <FC> <input.vcf> <output.vcf>

HELP
	exit;
}

# parameters
my $cutoff = shift @ARGV;
my $file = shift @ARGV;
my $outfile = shift @ARGV;

# input file
my $Stream = Bio::ToolBox::Data::Stream->new(
	file      => $file,
) or die "unable to open file '$file'!";
die "not a vcf file!" unless $Stream->vcf;

# output file
unless ($outfile) {
	$outfile = $Stream->path . $Stream->basename . '.filter' . $Stream->extension;
}
my $Out = $Stream->duplicate($outfile);

# filter
my $passCount = 0;
my $failCount = 0;
while (my $row = $Stream->next_row) {
	my $attrib = $row->vcf_attributes;
	my $fc = $attrib->{INFO}{FOLD_CHANGE_LOG} || 0;
	if ( abs($fc) >= $cutoff) {
		$Out->write_row($row);
		$passCount++;
	}
	else {
		$failCount++;
	}
}

# report
print "passed $passCount variants, failed $failCount variants, wrote $outfile\n";
