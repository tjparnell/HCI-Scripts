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
use IO::File;
use File::Basename qw(fileparse);

######## BarCode search pattern ########
# change this as necessary
# I could implement as options but I'm too lazy right now

# the length of the bar code
my $randomLength = 5;

# the fixed portion of the bar code between the random and the genomic read sequence
my $fixed = 'CTGA';

########################################

unless (@ARGV) {
	print <<END;

A script to pull out ChIP-Nexus barcodes out of a fastq and store it in 
the read name. This is to facilitate de-duplication using the random 
barcode to avoid PCR-duplicates. 
The random portion of the bar code is appended to the read name as ':NNNNN'.

Assumptions, limitations, and processing:
 - Bar codes are not checked for quality. 
 - Bar codes are not checked for mismatches, in either the random or fixed portion.
 - Reads that don't contain the bar code are skipped.
 - The number of processed and skipped are reported to stdout.
 - Assumes single-end fastq files. Paired-end files will get out of sync.
 - Assumes gzipped files. 
 - Takes as many input files as you give.
 - Output file is the input file basename appended with '.bc'. 

See the related script ChIPNexus_bam_dedup.pl for de-duplication.

Usage: $0 <file1.fastq.gz> <...>
END
	exit;
}

# calculate barcode length
my $barCodeLength = $randomLength + length($fixed);

# loop through the files
while (@ARGV) {
	my $file = shift @ARGV;
	print "processing $file....\n";
	my $goodCount = 0;
	my $badCount  = 0;
	
	# open input file
	my $infh = IO::File->new("gzip -dc $file|") or 
		die "unable to open $file! $!";
	
	# open output file
	my ($basename, $path, $extension) = 
		fileparse($file, qw(.txt.gz .fastq.gz));
	my $outfile = $path . $basename . '.bc' . $extension;
	my $outfh = IO::File->new("| gzip >$outfile") or 
		die "unable to open $outfile! $!";
	
	while (my $header  = $infh->getline) {
		chomp $header; 
		# we only chomp the header, we can safely skip the others
		my $sequence = $infh->getline or die "malformed file! no sequence line";
		my $spacer = $infh->getline or die "malformed file! no spacer line";
		my $quality = $infh->getline or die "malformed file! no quality line";
		
		# check for bar code sequence
		# I tried to set this up as a qr variable, but then nothing was being captured!?
		# so do it the old fashioned way
		if ($sequence =~ /^( [GATC] {$randomLength} ) $fixed/x) {
			my $random = $1 or die "nothing captured!"; # grabbed from the barcode regex
			my ($name, $desc) = split /\s+/, $header;
			$name .= ":$random";
			$header = $name . ' ' . $desc;
			$sequence = substr($sequence, $barCodeLength);
			$quality = substr($quality, $barCodeLength);
			$outfh->print("$header\n$sequence$spacer$quality");
			$goodCount++;
		}
		else {
			$badCount++;
		}
	}
	
	# close and finish
	$infh->close;
	$outfh->close;
	printf " %12s reads had a bar code and were processed\n %12s reads had no bar code\n", 
		$goodCount, $badCount;
	print " wrote $outfile\n";
}






