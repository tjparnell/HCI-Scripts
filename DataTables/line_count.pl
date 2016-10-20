#!/usr/bin/perl

use strict;
use Bio::ToolBox::Data;
use Bio::ToolBox::utility;

unless (@ARGV) {
	print <<HELP;

A quick line count program to count only the number of data lines, 
skipping comment, metadata, and blank lines. Handles gzipped files, 
BED, GFF, UCSC, VCF, and other text files. More than 2 input files 
are counted in parallel; output is not sorted.

Usage: $0 <file1> <file2> ...
HELP
	exit;
}

if (scalar @ARGV >= 2) {
	require Parallel::ForkManager;
	my $pm = Parallel::ForkManager->new(4);
	foreach my $f (@ARGV) {
		$pm->start and next;
		my $Data = Bio::ToolBox::Data->new(file => $f);
		printf "%10s %s\n", format_with_commas($Data->last_row), $f;
		$pm->finish;
	}
	$pm->wait_all_children;
	exit;
}
else {
	foreach my $f (@ARGV) {
		my $Data = Bio::ToolBox::Data->new(file => $f);
		printf "%12s %s\n", format_with_commas($Data->last_row), $f;
	}
}
