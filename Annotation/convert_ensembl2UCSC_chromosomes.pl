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
my $VERSION = 1.3;

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

Supports any recognizable annotation table (gtf gff3 bed refFlat genePred vcf)
and fasta files.

Output file is optional; it will append .ucsc to the input basename.

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
if ($infile =~ /\.fa(?:sta)?(?:\.gz)?$/i) {
	# input is a fasta file
	process_fasta();
}
else {
	# otherwise assume some sort of gene table
	process_table();
}




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

sub process_table {
	my $Stream = Bio::ToolBox::Data->new(
		stream => 1,
		in     => $infile,
	) or die " unable to open input table $infile! $!\n";
	# we check the chromosome column below
	
	# make output
	unless ($outfile) {
		$outfile = $Stream->path . $Stream->basename . '.ucsc' . $Stream->extension;
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
				my $chr2 = $chr;
				$chr2 =~ s/\.\d$//;
				$chr2 =~ s/\.\d$//;
				if (exists $lookup->{lc $chr2}) {
					my $alt = $lookup->{lc $chr2};
					$c =~ s/$chr/$alt/;
				}
			}
			elsif ($c =~ /^##contig=<ID=([\w\.]+)/) {
				# vcf sequence identifiers
				my $chr = $1;
				my $chr2 = $chr;
				$chr2 =~ s/\.\d$//;
				if (exists $lookup->{lc $chr2}) {
					my $alt = $lookup->{lc $chr2};
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

sub process_fasta {
	# open generic file handle to the input fasta file
	my $fh = Bio::ToolBox::Data->open_to_read_fh($infile) or 
		die " unable to open input file '$infile' $! \n";
	
	# open output file handle
	unless ($outfile) {
		$outfile = $infile;
		$outfile =~ s/\.fa(sta)?/.ucsc.fasta/g;
	}
	my $out = Bio::ToolBox::Data->open_to_write_fh($outfile) or
		die " unable to open output file '$outfile' $!\n";
	
	# conversion
	my %notfound;
	while (my $line = $fh->getline) {
		if (substr($line,0,1) eq '>') {
			# a fasta header line
			chomp($line); # because regex usually removes it anyway
			if ($line =~ /^>([\w\-\.]+)(\s+.+)?$/) {
				my $chr = $1;
				my $desc = $2;
				if (exists $lookup->{lc $chr}) {
					$chr = $lookup->{lc $chr};
				}
				else {
					$notfound{$chr}++;
				}
				$out->print(">$chr$desc\n");
			}
			else {
				warn "skipping malformed fasta header line '$line'\n";
				$out->print("$line\n");
			}
		}
		else {
			# a sequence line
			$out->print($line);
		}
	}
	$fh->close;
	$out->close;
	printf "wrote %s\n", $outfile;
	if (%notfound) {
		printf "could not convert the following chromosomes:\n%s\n", 
			join("\n", keys %notfound);
	}
}




