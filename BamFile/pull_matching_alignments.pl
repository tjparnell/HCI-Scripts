#!/usr/bin/perl

use strict;
use IO::File;
use Bio::DB::Sam;


unless (@ARGV) {
	print <<USAGE;

Extract sam alignments based on name.

Usage: $0 <input_bam> <list> <output_bam>

  <input_bam> is the source bam file
  <list> is a text file containing a list of read (query) names, one name per line.
  <output_bam> is the output bam file. This is optional

Alignments from the input bam file with a matching name are written to a new bam file.
USAGE
	exit 0;
}

# arguments
my $infile   = shift @ARGV;
my $namefile = shift @ARGV;
my $outfile  = shift @ARGV;

# check for required
die "need a list file!\n" unless $namefile;
unless ($outfile) {
	$outfile = $infile;
	$outfile =~ s/\.bam$/_matched.bam/i;
}

# load up list file
my %lookup;
my $duplicates = 0;
my $fh = IO::File->new($namefile);
while (my $line = $fh->getline) {
	chomp $line;
	if (exists $lookup{$line}) {
		$duplicates++;
	}
	else {
		$lookup{$line} = 1;
	}
}
$fh->close;
print " Warning: there were $duplicates duplicate query names in the list file!\n" if $duplicates;
printf "Loaded %s query names\n", scalar keys %lookup;

# open bam files
my $in = Bio::DB::Bam->open($infile) or die " Cannot open input Bam file!\n";
my $out = Bio::DB::Bam->open($outfile, 'w') or 
	die "cannot open output bam file '$outfile'!\n";

# write header
$out->header_write( $in->header );

# look for alignments
my $matched = 0;
while (my $a = $in->read1) {
	if (exists $lookup{ $a->qname }) {
		$matched++;
		$out->write1($a);
	}
}

# check for leftovers
if (%lookup) {
	printf "Warning: %s names could not be found!\n", scalar keys %lookup;
}

print "Matched $matched alignments\n";
print "Wrote output file $outfile\n";
exit 0;



