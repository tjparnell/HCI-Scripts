#!/usr/bin/perl

use strict;
use Bio::ToolBox::db_helper::bam; 

unless (@ARGV) {
	print <<END;

A script to remove excessive duplicate alignments, leaving up to X number 
of alignments at any give position, where X can be assigned by the user. 
This is in contrast to traditional duplicate removers that retain only 
one alignment per position.

Currently works with single-end data, for now. Paired-end alignments are 
treated like single end, likely breaking pairs.

Alignments are not selected for retention; the first X are retained.

Usage: $0 <limit> <input.bam> <output.bam>

END
	exit;
}

# max limit of reads at any given position
my $limit = shift @ARGV;
die "$limit is not an integer!\n" unless $limit =~ /^\d+$/;

# input bam file
my $infile = shift @ARGV;
my $sam = open_bam_db($infile) or # automatically takes care of indexing
	die "unable to read bam $infile $!";

# output bam file
my $outfile = shift @ARGV or die "no output file provided!\n";
my $outbam = write_new_bam_file($outfile) or 
	die "unable to open output bam file $outfile! $!";
	# this uses low level Bio::DB::Bam object

# write header
my $header = $sam->bam->header;
$outbam->header_write($header);

#initialize counters
my $totalCount = 0;
my $posCount   = 0;
my $keepCount  = 0;
my $tossCount  = 0;
my %dup2count; # for generating a histogram of duplicate numbers
               # key=duplicate count, value=number of bases with this rate

# process on chromosomes
# this could be parallelized per chromosome to speed up, if this is too slow
for my $tid (0 .. $sam->n_targets - 1) {
	my $seq_length = $sam->target_len($tid);
	
	# prepare callback data structure
	# anonymous hash of 
	my $data = {
		position   => 0,
		reads      => [],
	};
	
	# walk through the reads on the chromosome
	$sam->bam_index->fetch($sam->bam, $tid, 0, $seq_length, \&callback, $data);
	
	# check to make sure we don't leave something behind
	if ($data->{reads}->[0]) {
		write_reads($data);
	}
}

# print results
printf "
 %12s total mapped reads
 %12s mapped reads retained
 %12s duplicate mapped reads discarded
 original duplicate rate %.3f
 new duplicate rate %.3f
", $totalCount, $keepCount, $tossCount, ($totalCount - $posCount)/$totalCount, 
	$tossCount/$totalCount;
# print " ReadNumber NumberPositions\n";
# foreach (sort {$a <=> $b} keys %dup2count) {
# 	printf "%4s        %d\n", $_, $dup2count{$_};
# }

# finish up
undef $outbam;
check_bam_index($outfile);


### alignment callback
sub callback {
	my ($a, $data) = @_;
	$totalCount += 1;
	
	# check position 
	my $pos = $a->pos;
	if ($pos == $data->{position}) {
		push @{ $data->{reads} }, $a;
	}
	else {
		write_reads($data);
		$data->{position} = $pos;
		push @{ $data->{reads} }, $a;
	}
}

### write passing alignments 
sub write_reads {
	my $data = shift;
	
	# record how many at this position
	my $number = scalar @{ $data->{reads} };
	$dup2count{$number} += 1;
	$posCount++;
	
	# write out the max number of reads
	for (my $i = 0; $i < $limit; $i++) {
		my $a = shift @{ $data->{reads} };
		last unless $a;
		$outbam->write1($a);
		$keepCount++;
	}
	
	# record how many were lost
	$tossCount += scalar @{ $data->{reads} };
	
	# reset
	$data->{reads} = [];
}


