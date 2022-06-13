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
# https://github.com/tjparnell/HCI-Scripts


use strict;
use IO::File;
use List::Util qw(sum0);

my $VERSION = 1.0;

unless (@ARGV) {
	print <<END;
A script to combine Novoalign, Samtools markedup, and USeq bisulfite reports.
It scans all text files (extension .txt only) in each provided directory and 
extracts pertinant numbers, then combines them into a single output text file.
Samples are organized by column in the output file, suitable for import into 
Excel. Values not found are simply reported as zero.

Version: $VERSION

Usage: 
    
    combine_methyl_useq_reports.pl <outputfile> dir1/ dir2/ ....
    
END
	exit;
}

# output file
my $outfile = shift @ARGV;
if ($outfile !~ /\.txt$/i) {
	$outfile .= '.txt';
}

# process input files
my @dirs = @ARGV;
unless (@dirs) {
	die " no directories provided!\n";
}

my %samples;
foreach my $dir (@dirs) {
	unless (-d $dir) {
		print " $dir is not a directory! skipping\n";
		next;
	}
	print " scanning $dir\n";
	process_directory($dir);
}






# generate sort order
my @order;
{
	# tease apart the names
	my %sorter = (
		num  => {},
		char => {}
	);
	foreach my $path (@dirs) {
		if ($path =~ /^(\d+)X(\d+)/) {
			# GNomEx identifier
			$sorter{num}{$1}{$2} = $path;
		}
		else {
			$sorter{char}{$path} = $path;
		}
	}
	# sort by experiment ID, requestXsample, e.g. 1234X1
	foreach my $i (sort {$a <=> $b} keys %{$sorter{num}}) {
		foreach my $j (sort {$a <=> $b} keys %{$sorter{num}{$i}}) {
			# check to see if we have an array of multiple samples
			push @order, $sorter{num}{$i}{$j};
		}
	}
	# sort everything else
	foreach my $i (sort {$a cmp $b} keys %{$sorter{char}}) {
		push @order, $sorter{num}{$i};
	}
}


# output file
my $fh = IO::File->new($outfile, 'w') or die "can't write to $outfile!\n";

# header
$fh->printf("%s\n", join("\t", 'Program', 'Metric', @order));

# Novoalign metrics
$fh->printf("Novoalign\tTotal Reads\t%s\n", join("\t",
	map { $samples{$_}{novoalign_read} } @order));
$fh->printf("Novoalign\tUnique Mapped\t%s\n", join("\t",
	map { $samples{$_}{novoalign_unique} } @order));
$fh->printf("Novoalign\tUnique Mapped Percent\t%s\n", join("\t",
	map { $samples{$_}{novoalign_unique_percent} } @order));
$fh->printf("Novoalign\tMulti Mapped\t%s\n", join("\t",
	map { $samples{$_}{novoalign_multi} } @order));
$fh->printf("Novoalign\tUnique Mapped Percent\t%s\n", join("\t",
	map { $samples{$_}{novoalign_multi_percent} } @order));
$fh->printf("Novoalign\tUnique Unmapped\t%s\n", join("\t",
	map { $samples{$_}{novoalign_unmapped} } @order));
$fh->printf("Novoalign\tUnique Unmapped Percent\t%s\n", join("\t",
	map { $samples{$_}{novoalign_unmapped_percent} } @order));
$fh->printf("Novoalign\tPaired Insert Mean Size\t%s\n", join("\t",
	map { $samples{$_}{novoalign_insert_mean} } @order));

# Samtools markduplicates
$fh->printf("Samtools Markdup\tRead\t%s\n", join("\t",
	map { $samples{$_}{alignments_READ} } @order));
$fh->printf("Samtools Markdup\tWritten\t%s\n", join("\t",
	map { $samples{$_}{alignments_WRITTEN} } @order));
$fh->printf("Samtools Markdup\tExcluded\t%s\n", join("\t",
	map { $samples{$_}{alignments_EXCLUDED} } @order));
$fh->printf("Samtools Markdup\tPaired\t%s\n", join("\t",
	map { $samples{$_}{alignments_PAIRED} } @order));
$fh->printf("Samtools Markdup\tSingle\t%s\n", join("\t",
	map { $samples{$_}{alignments_SINGLE} } @order));
$fh->printf("Samtools Markdup\tCoordinate Primary Duplicates\t%s\n", join("\t",
	map { sum0($samples{$_}{alignments_DUP_PAIR}, $samples{$_}{alignments_DUP_SINGLE}) } @order));
$fh->printf("Samtools Markdup\tOptical Primary Duplicates\t%s\n", join("\t",
	map { sum0($samples{$_}{alignments_OPT_PAIR}, $samples{$_}{alignments_OPT_SINGLE}) } @order));
$fh->printf("Samtools Markdup\tTotal Duplicates\t%s\n", join("\t",
	map { $samples{$_}{alignments_DUP_TOTAL} } @order));
$fh->printf("Samtools Markdup\tDuplicate Fraction\t%s\n", join("\t",
	map {
		sprintf("%.3f", sum0(
			$samples{$_}{alignments_DUP_PAIR}, $samples{$_}{alignments_DUP_SINGLE},
			$samples{$_}{alignments_OPT_PAIR}, $samples{$_}{alignments_OPT_SINGLE}
		) / ($samples{$_}{alignments_READ} || 1))
	} @order) );

# USeq NovoalignBisulfiteParser 
$fh->printf("USeq NovoalignBisulfiteParser\tTotal alignments\t%s\n", join("\t",
	map { $samples{$_}{total_alignments} } @order));
$fh->printf("USeq NovoalignBisulfiteParser\tFailed Mapping quality\t%s\n", join("\t",
	map { $samples{$_}{alignments_failed_mapq} } @order));
$fh->printf("USeq NovoalignBisulfiteParser\tFailed Alignment score\t%s\n", join("\t",
	map { $samples{$_}{alignments_failed_AS} } @order));
$fh->printf("USeq NovoalignBisulfiteParser\tPassed\t%s\n", join("\t",
	map { $samples{$_}{alignments_passed} } @order));
$fh->printf("USeq NovoalignBisulfiteParser\tPassed Fraction\t%s\n", join("\t",
	map { $samples{$_}{alignments_passed_percent} } @order));
$fh->printf("USeq NovoalignBisulfiteParser\tNon-converted C Count\t%s\n", join("\t",
	map { $samples{$_}{nonconverted_C} } @order));
$fh->printf("USeq NovoalignBisulfiteParser\tConverted C count\t%s\n", join("\t",
	map { $samples{$_}{converted_C} } @order));
$fh->printf("USeq NovoalignBisulfiteParser\tFraction Non-converted Cs\t%s\n", join("\t",
	map { $samples{$_}{fraction_nonconverted} } @order));
$fh->printf("USeq NovoalignBisulfiteParser\tFraction BP passing quality\t%s\n", join("\t",
	map { $samples{$_}{fraction_bp_pass_quality} } @order));

# USeq BisStat 
$fh->printf("USeq BisStat\tLambda Conversion Rate\t%s\n", join("\t",
	map { $samples{$_}{labmda_fraction_nonc} } @order));
$fh->printf("USeq BisStat\tFraction Methyl C\t%s\n", join("\t",
	map { $samples{$_}{fraction_methyl_C} } @order));
$fh->printf("USeq BisStat\tFraction Methyl CpG\t%s\n", join("\t",
	map { $samples{$_}{fraction_methyl_CpG} } @order));

$fh->close;
printf "combined %d samples into $outfile\n", scalar(@order);



sub process_directory {
	
	my $dir = shift;
	$dir =~ s/\/$//; # remove any trailing slash
	
	my @files = glob "$dir/*.txt";
	
	# counts
	$samples{$dir} = {
		total_alignments            => 0,
		alignments_failed_mapq      => 0,
		alignments_failed_AS        => 0,
		alignments_passed           => 0,
		alignments_passed_percent   => 0,
		nonconverted_C              => 0,
		converted_C                 => 0,
		fraction_nonconverted       => 0,
		fraction_bp_pass_quality    => 0,
		labmda_fraction_nonc        => 0,
		fraction_methyl_C           => 0,
		fraction_methyl_CpG         => 0,
		alignments_READ             => 0,
		alignments_WRITTEN          => 0,
		alignments_EXCLUDED         => 0,
		alignments_PAIRED           => 0,
		alignments_SINGLE           => 0,
		alignments_DUP_PAIR         => 0,
		alignments_DUP_SINGLE       => 0,
		alignments_OPT_PAIR         => 0,
		alignments_OPT_SINGLE       => 0,
		alignments_DUP_TOTAL        => 0,
		novoalign_read              => 0,
		novoalign_unique            => 0,
		novoalign_unique_percent    => 0,
		novoalign_multi             => 0,
		novoalign_multi_percent     => 0,
		novoalign_unmapped          => 0,
		novoalign_unmapped_percent  => 0,
		novoalign_insert_mean       => 0
	};


	
	# process files
	foreach my $file (@files) {
		print "   scanning $dir/$file....\n"; 
		my $fh = IO::File->new($file) or 
			die "unable to read $file! $!\n";
	
	
		while (my $line = $fh->getline) {
			
			# USeq NovoalignBisulfiteParser
			if ($line =~  /^Filtering statistics for (\d+) alignments:$/) {
				$samples{$dir}{total_alignments} = $1;
			}
			elsif ($line =~ /^(\d+)\tFailed mapping quality score /) {
				$samples{$dir}{alignments_failed_mapq} = $1;
			}
			elsif ($line =~ /^(\d+)\tFailed alignment score /) {
				$samples{$dir}{alignments_failed_AS} = $1;
			}
			elsif ($line =~ /^(\d+)	Passed filters \((\d\d\.?\d?)%\)$/) {
				$samples{$dir}{alignments_passed} = $1;
				$samples{$dir}{alignments_passed_percent} = $2;
			}
			elsif ($line =~ /^(\d+)\tTotal non\-converted Cs sequenced$/) {
				$samples{$dir}{nonconverted_C} = $1;
			}
			elsif ($line =~ /^(\d+)\tTotal converted Cs sequenced$/) {
				$samples{$dir}{converted_C} = $1;
			}
			elsif ($line =~ /^(\d\.\d+)\tFraction non converted C's\.$/) {
				$samples{$dir}{fraction_nonconverted} = $1;
			}
			elsif ($line =~ /^(\d\.\d+)\tFraction bp passing quality \(\d+\)$/) {
				$samples{$dir}{fraction_bp_pass_quality} = $1;
			}
			
			# USeq BisStat 
			elsif ($line =~ /^Using Lambda data to set the expected fraction non-converted Cs to (0\.\d+) /) {
				$samples{$dir}{labmda_fraction_nonc} = $1;
			}
			elsif ($line =~ /^\t(\d\.\d\d\d)\t\(\d+\/\d+\)\tmC\/\(C\+mC\)$/) {
				$samples{$dir}{fraction_methyl_C} = $1;
			}
			elsif ($line =~ /^\t(\d\.\d\d\d)\t\(\d+\/\d+\)\tmCG\/\(CG\+mCG\)$/) {
				$samples{$dir}{fraction_methyl_CpG} = $1;
			}
			
			# Samtools MarkDup
			elsif ($line =~ /^READ: (\d+)$/) {
				$samples{$dir}{alignments_READ} = $1;
			}
			elsif ($line =~ /^WRITTEN: (\d+)$/) {
				$samples{$dir}{alignments_WRITTEN} = $1;
			}
			elsif ($line =~ /^EXCLUDED: (\d+)$/) {
				$samples{$dir}{alignments_EXCLUDED} = $1;
			}
			elsif ($line =~ /^PAIRED: (\d+)$/) {
				$samples{$dir}{alignments_PAIRED} = $1;
			}			
			elsif ($line =~ /^SINGLE: (\d+)$/) {
				$samples{$dir}{alignments_SINGLE} = $1;
			}
			elsif ($line =~ /^DUPLICATE PAIR: (\d+)$/) {
				$samples{$dir}{alignments_DUP_PAIR} = $1;
			}
			elsif ($line =~ /^DUPLICATE SINGLE: (\d+)$/) {
				$samples{$dir}{alignments_DUP_SINGLE} = $1;
			}
			elsif ($line =~ /^DUPLICATE PAIR OPTICAL: (\d+)$/) {
				$samples{$dir}{alignments_OPT_PAIR} = $1;
			}
			elsif ($line =~ /^DUPLICATE SINGLE OPTICAL: (\d+)$/) {
				$samples{$dir}{alignments_OPT_SINGLE} = $1;
			}
			elsif ($line =~ /^DUPLICATE TOTAL: (\d+)$/) {
				$samples{$dir}{alignments_DUP_TOTAL} = $1;
			}
			
			# Novoalign standard error output
			elsif ($line =~  /^#\s+Read Sequences:\s+(\d+)/) {
				$samples{$dir}{novoalign_read} = $1;
			}
			elsif ($line =~  /^#\s+Unique Alignment:\s+(\d+) \( ?(\d+\.\d)%\)$/) {
				$samples{$dir}{novoalign_unique} = $1;
				$samples{$dir}{novoalign_unique_percent} = $2;
			}
			elsif ($line =~ /^#\s+Multi Mapped:\s+(\d+) \( ?(\d+\.\d)%\)$/) {
				$samples{$dir}{novoalign_multi} = $1;
				$samples{$dir}{novoalign_multi_percent} = $2;
			}
			elsif ($line =~ /^#\s+No Mapping Found:\s+(\d+) \( ?(\d+\.\d)%\)$/) {
				$samples{$dir}{novoalign_unmapped} = $1;
				$samples{$dir}{novoalign_unmapped_percent} = $2;
			}
			elsif ($line =~ /^#\s+Mean\s+(\d+),\s+Std Dev\s+(\d+\.\d)$/) {
				$samples{$dir}{novoalign_insert_mean} = $1;
			}
			
		}
		$fh->close;
	}

}


