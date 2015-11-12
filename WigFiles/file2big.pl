#!/usr/bin/perl

# convert files to UCSC big formats

use strict;
use Bio::ToolBox::db_helper qw(open_db_connection);
use Bio::ToolBox::big_helper qw(wig_to_bigwig_conversion bed_to_bigbed_conversion);

my $help = <<HELP;
Convert files to equivalent UCSC big file formats
  Usage $0 <database> <file1> <file2> ...
  supports wig, bdg, bed
  database can be any indexed file (bam,bw,bb), fasta, or Bio::DB::SeqFeature::Store db
HELP

unless (@ARGV) {
	print $help;
	exit;
}

my $database = shift @ARGV;
my $db = open_db_connection($database) or 
	die "unable to open database '$database'!\n";

while (@ARGV) {
	my $file = shift @ARGV;
	my $b;
	if ($file =~ /(?:wig|bdg|bedgraph)$/i) {
		$b = wig_to_bigwig_conversion( db => $db, wig => $file );
	}
	elsif ($file =~ /bed$/) {
		$b = bed_to_bigbed_conversion( db => $db, bed => $file );
	}
	elsif ($file =~ /(?:gz|bz2|zip)$/) {
		print "please decompress $file first\n";
		next;
	}
	else {
		print "unrecognized extension for '$file'\n";
		next;
	}
	
	if ($b) {
		print "successfully converted $file to $b\n";
	}
	else {
		print "failed to convert $file\n";
	}
}



