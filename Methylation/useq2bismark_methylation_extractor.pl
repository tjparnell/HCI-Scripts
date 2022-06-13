#!/usr/bin/env perl

# a script to convert USeq CpG files into a Bismark compatible file format

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
use Bio::ToolBox::Data;
# must also have Bio::DB::USeq version 0.24 installed

unless (@ARGV) {
print <<END;
  
  Extract BisMark-like coverage file from ConvertedC and NonConvertedC USeq files
  from USeq NovoalignBisulfiteParser output files.
  
  Be sure to first filter point data for CpG context using the ParsePointDataContexts.
  It will silently skip mitochondrial and Lambda chromosomes.
  
  It will write .bismark.cov and .bedGraph files.
  
  Usage: bismark_methylation_extractor.pl <ConvertedCpG.useq> <NonConvertedCpG.useq> <out>

END
	exit;
}

my ($c_file, $nonc_file, $outfile) = @ARGV;

my $c_db = Bio::ToolBox::Data->open_new_database($c_file) or 
	die "cannot open $c_file! $!\n";
my $nonc_db = Bio::ToolBox::Data->open_new_database($nonc_file) or 
	die "cannot open $nonc_file! $!\n";
my $outfh1 = Bio::ToolBox::Data->open_to_write_fh($outfile . '.bismark.cov') or 
	die "cannot open file $outfile! $!\n";
my $outfh2 = Bio::ToolBox::Data->open_to_write_fh($outfile . '.bedGraph') or 
	die "cannot open file $outfile! $!\n";

foreach my $chr ($c_db->seq_ids) {
	print " processing $chr\n";
	
	# each position will be an array of [c, nonc]
	my %pos2score;
	
	# try to do this in as efficient manner as possible
	# collect the segement for this chromosome
	# identify the slices for this segment
	# collect the observations from each slice... 
	# Easier than collecting observations across the entire chromosome at once
	# and avoids wrapping each observation into a SeqFeature object that the 
	# feature iterator will do
	
	# need to process each strand separately
	# the output Bismark file format combines the C values from each strand 
	# for a CpG pair into a single value
	
	
	#### first collect converted Cs
	my $segment = $c_db->segment(-seq_id => $chr);
	foreach my $slice ($segment->slices) {
		# plus strand first
		my $obs = $c_db->observations(
			-seq_id   => $c_db->slice_seq_id($slice),
			-start    => $c_db->slice_start($slice),
			-end      => $c_db->slice_end($slice),
			-strand   => 1
		);
		foreach my $o (@$obs) {
			# each observation is array ref of start0, stop, score
			my $p = $o->[1];
			$pos2score{$p} ||= [0,0];
			$pos2score{$p}->[0] += $o->[2];
		}
		# minus strand next
		my $obs = $c_db->observations(
			-seq_id   => $c_db->slice_seq_id($slice),
			-start    => $c_db->slice_start($slice),
			-end      => $c_db->slice_end($slice),
			-strand   => -1 
		);
		foreach my $o (@$obs) {
			# each observation is array ref of start0, stop, score
			my $p = $o->[0];
			$pos2score{$p} ||= [0,0];
			$pos2score{$p}->[0] += $o->[2];
		}
	}
	
	#### Collect Nonconverted Cs
	$segment = $nonc_db->segment(-seq_id => $chr);
	foreach my $slice ($segment->slices) {
		# plus strand first
		my $obs = $nonc_db->observations(
			-seq_id   => $nonc_db->slice_seq_id($slice),
			-start    => $nonc_db->slice_start($slice),
			-end      => $nonc_db->slice_end($slice),
			-strand   => 1
		);
		foreach my $o (@$obs) {
			# each observation is array ref of start0, stop, score
			my $p = $o->[1];
			$pos2score{$p} ||= [0,0];
			$pos2score{$p}->[1] += $o->[2];
		}
		# minus strand next
		my $obs = $nonc_db->observations(
			-seq_id   => $nonc_db->slice_seq_id($slice),
			-start    => $nonc_db->slice_start($slice),
			-end      => $nonc_db->slice_end($slice),
			-strand   => -1
		);
		foreach my $o (@$obs) {
			# each observation is array ref of start0, stop, score
			my $p = $o->[0];
			$pos2score{$p} ||= [0,0];
			$pos2score{$p}->[1] += $o->[2];
		}
	}
	
	
	### write out
	foreach my $pos (sort {$a <=> $b} keys %pos2score) {
		my ($c, $nc) = @{ $pos2score{$pos} };
		my $frac = ($nc / ($c + $nc)) * 100;
		# coverage file
		$outfh1->printf("$chr\t$pos\t$pos\t%.0f\t%d\t%d\n", $frac, $nc, $c);
		# bedgraph file
		$outfh2->printf("$chr\t%d\t$pos\t%.0f\n", $pos - 1, $frac);
	}
}

$outfh1->close;
$outfh2->close;

print " finished\n";
