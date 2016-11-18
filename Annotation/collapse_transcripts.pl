#!/usr/bin/perl

use strict;
use Bio::ToolBox::Data;
use Bio::ToolBox::GeneTools 1.50 qw(collapse_transcripts ucsc_string);

unless (@ARGV) {
	print <<HELP;
A script to calculate collapse transcripts.

Will take any gene table format as input: gff3, gtf, ucsc tables.

Collapses genes into a single model and outputs a simple refFlat file.

Usage: $0 <table.gff> <collapsed.refFlat>

HELP
	exit;
}

my $start = time;

## Open input file
my $file = shift @ARGV;
my $Data = Bio::ToolBox::Data->new(file => $file, parse => 1) or 
	die "unable to open input file '$file'\n";
printf "Parsed %d genes in %.02f minutes\n", $Data->last_row, (time - $start)/60;

## Open output file
my $outfile = shift @ARGV or die "no output file given!\n";
my $outfh = $Data->open_to_write_fh($outfile);
$outfh->print("# Collapsed from $file with $0\n");

## Collapse and write
my $count = 0;
$Data->iterate( sub {
	my $row = shift;
	my $gene = $row->seqfeature;
	my $tiny = collapse_transcripts($gene);
	my $ucsc = ucsc_string($tiny);
	$outfh->print("$ucsc\n");
	$count++;
} );

$outfh->close;
printf "Collapsed $file into $count minimal gene models in %.02f minutes\n", (time - $start)/60;


