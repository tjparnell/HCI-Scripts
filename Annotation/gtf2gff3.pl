#!/usr/bin/perl

# convert GTF to GFF3

use strict;
use Bio::ToolBox::parser::gff;
use Bio::ToolBox::Data::file;

my $help = <<HELP;
Convert GTF files to GFF3 format
  Usage $0 <file1.gtf> <file2.gtf.gz> ...
HELP

unless (@ARGV) {
	print $help;
	exit;
}

while (@ARGV) {
	# process file
	my $file = shift @ARGV;
	print " parsing $file....\n";
	my $parser = Bio::ToolBox::parser::gff->new($file);
	my $p = $parser->parse_file;
	unless ($p) {
		"   failed!\n";
		next;
	}
	
	# prepare output
	my $out = $file;
	$out =~ s/gtf/gff3/i;
	my $fh = Bio::ToolBox::Data::file->open_to_write_fh($out); # preserves compression
	unless ($fh) {
		"  unable to write $out!\n";
		next;
	}
	
	# print metadata
	$fh->print("##gff-version 3\n");
	$fh->print("# converted from $file\n");
	foreach ($parser->comments) {
		$fh->print("$_\n");
	}
	
	# print features
	while (my $f = $parser->next_top_feature) {
		$f->version(3);
		$fh->print( $f->gff_string(1) ); # recursive gff printing
		$fh->print("\n###\n"); # close subfeature pragma
	}
	print "  wrote file $out\n";
}

