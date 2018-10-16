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
use Getopt::Long qw(:config no_ignore_case);
use List::Util qw(sum);
use Bio::ToolBox::Data;
use Bio::ToolBox::GeneTools 1.63 qw(
	:export
	:cds 
	:transcript
);
use Bio::ToolBox::utility;
my $VERSION = 1;

my $DOC = <<DOC;
  
  A simple script to trim inwards translation start and stop coordinates 
  of an annotation file, such as a GTF file. Note that CDS information must 
  be present in the file. Non-coding genes won't be affected. The effect 
  will be that CDS length is shorter. This may be useful for applications 
  such as ribosome profiling.
  
  USAGE:
  trim_cds_ends.pl --in file.gtf --begin 9 --end 12 --out file2.gtf
  
    -i --in <file>          Input GTF, GFF3, refFlat annotation file
    -b --begin <integer>    Number of bp to subtract from translation start
    -e --end <integer>      Number of bp to subract from translation end
    -o --out <file>         Output GTF file
    
DOC

unless (@ARGV) {
	print $DOC;
	exit;
}


#### Options
my $input;
my $start;
my $stop;
my $outfile;

GetOptions( 
	'i|in=s'      => \$input, # input table
	'b|begin=i'   => \$start, # how much to trim from start
	'e|end=i'     => \$stop, # how much to trim from end
	'o|out=s'     => \$outfile, # output file
) or die " unrecognized options!\n";

# minimum length 
my $min_length = ($start || 0) + ($stop || 0) + 3; # make it a little longer than requested
die " trim values add up to zero!!??\n" if $min_length == 0;
print " minimum length of CDS must be $min_length bp\n";


#### Input file
# we're using a Data object for simplicity
# this has a side benefit of working with multiple annotation file types
my $Data = Bio::ToolBox::Data->new(
	file       => $input, 
	parse      => 1,
	simplify   => 0, # we want everything!
	feature    => 'gene', # take all top level genes
	subfeature => 'exon,cds,utr',
) or die " unable to load input file '$input'\n";

if ($Data->last_row) {
	printf " Loaded %s genes from $input.\n", format_with_commas($Data->last_row);
}
else {
	die " No features loaded!\n";
}


#### Output
# open a simple output file handle
# limiting to GTF here only, complain if you need more
my $outfh = Bio::ToolBox::Data->open_to_write_fh($outfile) or 
	die "unable to open '$outfile' for writing! $!\n";
$outfh->print("##gff-version 2.5\n");
$outfh->print("# translation start trimmed by $start bp\n");
$outfh->print("# translation end trimmed by $stop bp\n");
$outfh->printf("# exported from %s\n", $input);



#### Iterate
my $count = 0;
$Data->iterate(\&callback);
$outfh->close;
print " wrote $outfile\n";
printf " Trimmed %s transcripts\n", format_with_commas($count);


# primary Data callback routine for working on each Data row, i.e. gene
sub callback {
	my $row = shift;
	my $gene = $row->seqfeature(1); # this should already be in memory
	
	# walk through each transcript
	foreach my $transcript (get_transcripts($gene)) {
		# skip non-coding
		next unless is_coding($transcript);
		
		# process based on strand
		if ($transcript->strand > 0) {
			adjust_forward_transcripts($transcript);
		}
		else {
			adjust_reverse_transcripts($transcript);
		}
		
		# remove UTRs and codons, since these will no longer be accurate
		# and I don't want to bother recalculating them right now
		# I'm not using GeneTools functions here, because they will automatically 
		# create these subfeatures if they're not present.
		foreach my $s ($transcript->get_SeqFeatures) {
			my $id = $s->primary_id; # force an ID to be generated if not present
			if ($s->primary_tag =~ /codon|utr/i) {
				$transcript->delete_SeqFeature($id);
			}
		} 
		
		$count++;
	}
	
	# now write the gene
	my $string = gtf_string($gene);
		# the gtf string will automatically recreate codons, but does not include UTRs
	$outfh->print($string);
}


sub adjust_forward_transcripts {
	my $transcript = shift;
	
	# get CDS
	my @cds = get_cds($transcript);
	my $cds_length = sum( map {$_->length} @cds );
	if ($cds_length <= $min_length) {
		printf " transcript %s CDS is $cds_length bp and too short to trim\n", $transcript->primary_id;
		return;
	}
	
	## First CDS
	# make sure first CDS is long enough
	if ($cds[0]->length > $start) {
		# it's long enough, so adjust
		my $s = $cds[0]->start + $start;
		$cds[0]->start($s);
	}
	else {
		my $remainder = $start - $cds[0]->length;
		$transcript->delete_SeqFeature( $cds[0]->primary_id );
		# adjust second, now first, CDS
		my $s = $cds[1]->start + $remainder;
		$cds[1]->start($s);
	}
	
	## Last CDS
	# make sure last CDS is long enough
	if ($cds[-1]->length > $stop) {
		# it's long enough, so adjust
		my $s = $cds[-1]->end - $stop;
		$cds[-1]->end($s);
	}
	else {
		my $remainder = $stop - $cds[-1]->length;
		$transcript->delete_SeqFeature( $cds[-1]->primary_id );
		# adjust second, now first, to last CDS
		my $s = $cds[-2]->end - $remainder;
		$cds[-2]->start($s);
	}
}


sub adjust_reverse_transcripts {
	my $transcript = shift;
	
	# get CDS
	my @cds = get_cds($transcript);
	my $cds_length = sum( map {$_->length} @cds );
	if ($cds_length <= $min_length) {
		printf " transcript %s CDS is $cds_length bp and too short to trim\n", $transcript->primary_id;
		return;
	}
	
	## First CDS
	# make sure first CDS is long enough
	if ($cds[-1]->length > $start) {
		# it's long enough, so adjust
		my $s = $cds[-1]->end - $start;
		$cds[-1]->end($s);
	}
	else {
		my $remainder = $start - $cds[-1]->length;
		$transcript->delete_SeqFeature( $cds[-1]->primary_id );
		# adjust second, now first, CDS
		my $s = $cds[-2]->end - $remainder;
		$cds[-2]->end($s);
	}
	
	## Last CDS
	# make sure last CDS is long enough
	if ($cds[0]->length > $stop) {
		# it's long enough, so adjust
		my $s = $cds[0]->start + $stop;
		$cds[0]->start($s);
	}
	else {
		my $remainder = $stop - $cds[0]->length;
		$transcript->delete_SeqFeature( $cds[0]->primary_id );
		# adjust second, now first, to last CDS
		my $s = $cds[1]->start + $remainder;
		$cds[1]->start($s);
	}
}

