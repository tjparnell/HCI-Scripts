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
# https://github.com/HuntsmanCancerInstitute/MultiRepMacsChIPSeq


use strict;
use File::Spec;
use Bio::ToolBox 1.68;
use Bio::ToolBox::db_helper 1.50 qw(
	open_db_connection 
	low_level_bam_fetch
);
use Bio::ToolBox::GeneTools qw(get_exons);

my $VERSION = 1;



unless (@ARGV) {
print <<USAGE;
  
  A script to identify cDNA alignments from a genomic bam file.
  In other words, contamination of an expression vector in otherwise 
  genomic DNA sequencing experiment. Checks just one gene against 
  one or more Bam files.
  
  Version: $VERSION
   
  Usage: $0 <output.txt> <gene_gtf> <file1.bam> <file2.bam ...>
  
USAGE
exit;
}

my $outfile = shift @ARGV;
my $annofile = shift @ARGV;
my @bams = @ARGV;


# Input annoation
my $Input = Bio::ToolBox->parse_file(
	file       => $annofile,
	simplify   => 1,
	feature    => 'gene',
	subfeature => 'exon,cds,utr'
);
if ($Input->last_row > 1) {
	warn " uh oh! your gtf file has more than one gene! just using the first one, ok?\n";
}
$Input->collapse_gene_transcripts; # just in case
my $stream = $Input->row_stream;
	# I just need the first row feature, but I have to make a stream object just to do it
my $feature = $stream->next_row; # a row feature 
my $gene = $feature->seqfeature;
my @exons = get_exons($gene);
unless (@exons) {
	die " uh oh! No exons for gene!\n";
}


# Output
my $Data = Bio::ToolBox->new_data( qw(File Mean_Exon_Coverage Mean_Intron_Coverage 
	Number_Genomic_Alignments Number_cDNA_Alignments));


# Go through bam files
foreach my $b (@bams) {
	
	my (undef, undef, $name) = File::Spec->splitpath($b);
	$name =~ s/\.bam$//i;
	print " processing $name...\n";
	
	# verify
	$b = $Input->verify_dataset($b);
	unless ($b) {
		warn "can't verify $b! skipping\n";
		next;
	}
	
	# classify alignments
	my $counter = {
		start   => 0,
		end     => 0,
		genomic => 0,
		cdna    => 0
	};
	my $sam = open_db_connection($b);
	my ($tid, undef, undef) = $sam->header->parse_region($feature->seq_id);
	foreach my $exon (@exons) {
		$counter->{start} = $exon->start - 1;
		$counter->{end} = $exon->end;
		low_level_bam_fetch($sam, $tid, $exon->start - 1, $exon->end, \&callback, $counter);
	}
	
	# exon coverage
	my $exon_cov = $feature->get_score(
		dataset    => $b,
		method     => 'mean',
		subfeature => 'exon',
	);
	
	# intron coverage
	my $intron_cov = $feature->get_score(
		dataset    => $b,
		method     => 'mean',
		subfeature => 'intron',
	);
	
	# record
	$Data->add_row( [
		$name,
		sprintf("%.2f", $exon_cov),
		sprintf("%.2f", $intron_cov),
		$counter->{genomic},
		$counter->{cdna}
	] );
	
}


# save
my $s = $Data->write_file($outfile);
print " wrote $s\n";


sub callback {
	my ($a, $data) = @_;
	my $cigar = $a->cigar_str;
	if ($cigar =~ /^\d+S/ and $a->pos == $data->{start}) {
		# left soft trim
		$data->{cdna} += 1;
	}
	elsif ($cigar =~ /\d+S$/ and $a->calend == $data->{end}) {
		# right soft trim
		$data->{cdna} += 1;
	}
	else {
		$data->{genomic} += 1;
	}
}






