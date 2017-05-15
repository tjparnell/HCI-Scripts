#!/usr/bin/env perl

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
use Getopt::Long;
use Bio::ToolBox::parser::ucsc '1.44';
use Bio::ToolBox::GeneTools 1.44 qw(gtf_string);
use Bio::ToolBox::utility qw(open_to_write_fh format_with_commas);

unless (scalar @ARGV) {
	print <<END;
A script to convert UCSC-style files to a GTF format.

Usage:
 $0 <infile> <outfile>
 
    infile may be any refFlat, genePred, or genePredExt format.
    outfile is optional, defaults to basename of input file.
END
	exit;
}


# input file
my $infile = shift (@ARGV);

# output file, optional
my $outfile = shift (@ARGV) || undef;
unless (defined $outfile) {
	$outfile = $infile;
	$outfile =~ s/(?:reff?lat|genepred|ucsc|txt)/gtf/i;
	print "assigning outfile to $outfile\n";
	die "oops! can't overwrite file with same name!\n" if $outfile eq $infile;
}

# open ucsc parser object
my $parser = Bio::ToolBox::parser::ucsc->new($infile) or
	die " unable to open input file '$infile'!\n";

# open output file handle
my $outfh = open_to_write_fh($outfile) or 
	die "unable to open output file '$outfile' $!\n";

# Process the top features
my @top_features = $parser->top_features();
printf " parsed %s top features\n", format_with_commas(scalar @top_features);
while (@top_features) {
	my $feature = shift @top_features;
	my $string = gtf_string($feature);
	$outfh->print($string);
}

$outfh->close;
print " wrote '$outfile'\n";






