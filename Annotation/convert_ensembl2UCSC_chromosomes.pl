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
use Bio::ToolBox::Data '1.41';
my $VERSION = 1.1;

unless (scalar @ARGV >= 2) {
	print <<END;

A script to convert Ensembl chromosome identifiers to UCSC identifiers.

Usage:
 $0 <chromInfo> <infile> <outfile>

chromInfo should be a list of the UCSC chromosomes. Download from UCSC.
Example:
  wget http://hgdownload.soe.ucsc.edu/goldenPath/mm10/database/chromInfo.txt.gz
  wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/chromInfo.txt.gz
  wget http://hgdownload.soe.ucsc.edu/goldenPath/rn6/database/chromInfo.txt.gz
You could also use a samtools fasta index .fai file. 

Supports any recognizable annotation table (gtf gff3 bed refFlat genePred vcf).

Output file is optional; it will append _ucsc to the input basename.

It will report any unmatched chromosomes that couldn't be converted. 
Check your chromosome file to see if these, or something like it, is 
really present. If so, contact the author for a fix.

END
	exit;
}

my $chromfile = shift @ARGV;

my $lookup = make_lookup_table();

my $infile = shift @ARGV;
my $outfile = shift @ARGV || undef;
process();




sub make_lookup_table {
	my $Data = Bio::ToolBox::Data->new() or die "no Data object!";
	my $fh = $Data->open_to_read_fh($chromfile) or 
		die "unable to open $chromfile! $!\n";
	my %lookup;
	while (my $line = $fh->getline) {
		chomp $line;
		my ($chr, undef) = split(/\s+/, $line, 2);
		my $alt;
		if ($chr =~ /^chr(\d+)$/) {
			$alt = lc $1;
		}
		elsif ($chr =~ /^chr[A-Za-z\d]+_([A-Za-z\d]+)v\d_?(?:random|alt)?$/) {
			$alt = lc $1;
		}
		elsif ($chr =~ /^chr[A-Za-z\d]+_([A-Za-z\d]+)_?(?:random|alt)?$/) {
			$alt = lc $1;
		}
		elsif ($chr eq 'chrM') {
			$alt = 'mt';
		}
		elsif ($chr =~ /^chr(\w+)$/) {
			$alt = lc $1;
		}
		else {
			$alt = lc $chr;
		}
		$lookup{$alt} = $chr;
# 		print " using $alt for $chr\n";
	}
	printf " using %d chromosomes in the lookup table\n", scalar keys %lookup;
	return \%lookup;
}

sub process {
	my $Stream = Bio::ToolBox::Data->new(
		stream => 1,
		in     => $infile,
	) or die " unable to open input table $infile! $!\n";
	# we check the chromosome column below
	
	# make output
	unless ($outfile) {
		$outfile = $Stream->path . $Stream->basename . '_ucsc' . $Stream->extension;
	}
	my $Out = $Stream->duplicate($outfile) or 
		die " unable to open output stream for file $outfile! $!\n";
	
	# deal with metadata
	my @comments = $Stream->comments;
	if (@comments) {
		for (my $i =$#comments; $i >= 0; $i--) {
			# delete the existing comments, these are indexed so go in reverse
			# order, we'll add back fixed ones
			$Out->delete_comment($i); 
		}
		foreach my $c (@comments) {
			if ($c =~ /^##sequence\-region\s+([\w\.]+)\s/) {
				# gff3 sequence pragmas
				my $chr = $1;
				$chr =~ s/\.\d$//;
				if (exists $lookup->{lc $chr}) {
					my $alt = $lookup->{lc $chr};
					$c =~ s/$chr/$alt/;
				}
			}
			elsif ($c =~ /^##contig=<ID=([\w\.]+)/) {
				# vcf sequence identifiers
				my $chr = $1;
				$chr =~ s/\.\d$//;
				if (exists $lookup->{lc $chr}) {
					my $alt = $lookup->{lc $chr};
					$c =~ s/$chr/$alt/;
				}
			}
			$Out->add_comment($c);
		}
		
	}
	
	# data replacements
	my $seq_i = $Stream->chromo_column;
	die "can't find chromosome column!\n" unless defined $seq_i;
	my %notfound;
	while (my $row = $Stream->next_row) {
		my $chr = $row->value($seq_i);
		$chr =~ s/\.\d$//;
		if (exists $lookup->{lc $chr}) {
			$row->value($seq_i, $lookup->{lc $chr});
		}
		else {
			$notfound{$chr}++;
		}
		$Out->write_row($row);
	}
	$Out->close_fh;
	$Stream->close_fh;
	printf "wrote %s\n", $Out->filename;
	if (%notfound) {
		printf "could not convert the following chromosomes:\n%s\n", 
			join("\n", keys %notfound);
	}
}
