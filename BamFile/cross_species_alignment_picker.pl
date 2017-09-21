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
use Getopt::Long;
use Bio::ToolBox::db_helper 1.50 qw(
	open_db_connection 
	$BAM_ADAPTER
);

my $VERSION = 1.3;
# version 1.0 used mapping quality score and number of cigar operations
# version 1.1 uses alignment score and sum of mismatches and cigar operations
# version 1.2 improves text output to make counts less confusing
# version 1.3 add support for other aligner AS tags and make bam adapter agnostic

my $description = <<DESCRIPTION;

Picks the best alignment between two species.

Sequence reads from a mixed species sample, e.g. human xenografts in mouse, 
are aligned to each species index. Unique species alignments are kept, while 
alignments matching to both species are either discarded, written to a cross 
species file, or assigned to one of the two species.

Provide two bam files representing alignments to each species index. Bam files 
do not need to be sorted, but any sort order will be maintained. Alignments 
marked as supplementary, secondary, or unmapped are silently discarded; 
everything else is retained. Paired-end alignments are not currently 
supported. 

Alignments mapping to both species may be picked as to the best alignment 
based on two tests of three criteria (one test isn't always enough). These  
tests include the following:
    1. Alignment score, stored as the AS attribute tag. This number is 
       aligner specific. In most aligners, the higher score is a better  
       alignment. The exception is Novoalign, where lower alignment 
       score is better.
    2. The number of mismatchs, stored as the NM attribute tag. This 
       is added to the number of CIGAR string operations, which can 
       include insertions, deletions, and trims. Presumably, a lower 
       number of mismatches and operations is a better alignment. 

Alignments equal in the metric score are considered orthologous.
Othologous reads mapping to both indices can either be discarded (default) 
or written to a second bam file. When the best pick option is enabled, 
the equivalently mapping files are written to the orthologous bam file.

For each input bam file, a .unique.bam file is written. If specified, 
orthologous files are written to a second output bam file, .orthologous.bam.

DESCRIPTION

my $usage = <<USAGE;

Usage: $0 [options] file1.bam file2.bam 

Options:
    --ortho     Write a second bam file with all orthologous or 
                equally mapping alignments
    --pick      Assign the best alignment to one species based on 
                the alignment test
    --low       The lower Alignment Score wins (Novoalign)

USAGE



### Options
	# ugh, GetOptions doesn't like it when no options are given
	# so gotta check if we have options or not
my $do_ortho = 0;
my $do_pick = 0;
my $low_score;
my ($file1, $file2);

if (scalar @ARGV == 0) {
	print $description . $usage;
	exit 0;
}
if (scalar @ARGV == 2) {
	($file1, $file2) = @ARGV;
}
elsif (scalar @ARGV > 2) {
	GetOptions( 
		'ortho!'     => \$do_ortho,
		'pick!'      => \$do_pick,
		'low!'       => \$low_score,
		'bam=s'      => \$BAM_ADAPTER,
	) or die "bad options!\n$usage";
	($file1, $file2) = @ARGV;
}
else {
	die "bad options!\n$usage";
}




### Global values
my %alignments;
my $count1best = 0;
my $count2best = 0;
my $countequal = 0;
my $countortho = 0;

# set the subroutine to collect the Alignment Score metric
my $get_metric1 = $low_score ? \&get_low_as_metric : \&get_high_as_metric;




### Check the alignments in both files
my $file1count = read_first_bam();
my $file2count = read_second_bam();

# report numbers
printf " There were %d mapped alignments in $file1\n", $file1count;
printf "   %d alignments mapped better than the second file\n", $count1best 
	if $do_pick;
printf " There were %d mapped alignments in $file2\n", $file2count;
printf "   %d alignments mapped better than the first file\n", $count2best 
	if $do_pick;
printf " There were %d equally mapped alignments\n", $countequal if $do_pick;

if ($do_ortho and not $do_pick) {
	printf " There were %d orthologous alignments\n", $countortho;
}




### Write out the reads
write_bam_files($file1, 1);
write_bam_files($file2, 2);

exit;




### Subroutines

sub read_first_bam {
	
	my $in = open_db_connection($file1) or die " Cannot open input Bam file!\n";
	print " Reading $file1...\n";
	my $count = 0;
	
	# header
	my $bam;
	my $header;
	if ($BAM_ADAPTER eq 'sam') {
		$bam = $in->bam;
		$header = $bam->header;
	}
	elsif ($BAM_ADAPTER eq 'hts') {
		$bam = $in->hts_file;
		$header = $bam->header_read;
	}
	else {
		die "unrecognized bam adapter $BAM_ADAPTER!";
	}
	
	# walk through reads
	while (my $a = $bam->read1($header)) {
		next if $a->unmapped;
		my $flag = $a->flag;
		next if ($flag & 0x0100); # secondary alignment
		next if ($flag & 0x0800); # supplementary hit
		my $score1 = $get_metric1->($a);
		my $score2 = get_mismatch_metric($a);
		$alignments{$a->qname} = "$score1,$score2";
		$count++;
	}
	return $count;
}

sub read_second_bam {
	# this reads the second bam file, and compares the metric 
	# number with every mapped alignment that matches from the first file
	# alignments are assigned a single value 
	# 1 is keep first file read
	# 2 is keep second file read
	# 3 is orthologous read in both, keep or discard as requested
	# 0 is discard (multiple hit alignment)
	
	my $in = open_db_connection($file2) or die " Cannot open input Bam file!\n";
	print " Reading $file2...\n";
	my $count = 0;

	# header
	my $bam;
	my $header;
	if ($BAM_ADAPTER eq 'sam') {
		$bam = $in->bam;
		$header = $bam->header;
	}
	elsif ($BAM_ADAPTER eq 'hts') {
		$bam = $in->hts_file;
		$header = $bam->header_read;
	}
	else {
		die "unrecognized bam adapter $BAM_ADAPTER!";
	}
	
	# walk through reads
	while (my $a = $bam->read1($header)) {
		next if $a->unmapped;
		my $flag = $a->flag;
		next if ($flag & 0x0100); # secondary alignment
		next if ($flag & 0x0800); # supplementary hit
		my $name = $a->qname;
		$count++;
	
		# check for existing
		if (exists $alignments{$name}) {
			next if substr($alignments{$name},0,1) eq '='; 
				# we already processed this so move one
		
			# compare
			if ($do_pick) {
				my ($score1as, $score1mm) = split ',', $alignments{$name};
				my $score2as = $get_metric1->($a);
				if ($score1as > $score2as) {
					# first alignment is better
					$alignments{$name} = '=1';
					$count1best++;
				}
				elsif ($score1as < $score2as) {
					# second alignment is better
					$alignments{$name} = '=2';
					$count2best++;
				}
				elsif ($score1as == $score2as) {
					# alignments are equal
					# compare second score
					my $score2mm = get_mismatch_metric($a);
					
					if ($score1mm < $score2mm) {
						# first alignment is better
						$alignments{$name} = '=1';
						$count1best++;
					}
					elsif ($score1mm > $score2mm) {
						# second alignment is better
						$alignments{$name} = '=2';
						$count2best++;
					}
					else {
						# alignments are truly equal
						if ($do_ortho) {
							# keep for third orthologous bam file
							$alignments{$name} = '=3';
						}
						else {
							# drop the alignment
							$alignments{$name} = '=0';
						}
						$countequal++;
					}
				}
				next;
			}
			elsif ($do_ortho) {
				# keep for third orthologous bam file
				$alignments{$name} = '=3';
				$countortho++;
			}
			else {
				# discard it 
				$alignments{$name} = '=0' if exists $alignments{$name};
			}
		}
		else {
			# great, unique to the second file
			$alignments{$name} = '=2';
		}
	}
	return $count;
}


sub write_bam_files {
	# writes the unique and orthologous bam files for each input file
	# walks through input file again, and writes alignment to appropriate
	# output file or skips and moves on
	
	my ($file, $number) = @_;
	my $count = 0;
	my $cross_count = 0;
	
	# input file
	my $in = open_db_connection($file, 1); # force to open new file handle 
	
	# header and bam writer
	my $bam;
	my $header;
	my $write_alignment;
	if ($BAM_ADAPTER eq 'sam') {
		$bam = $in->bam;
		$header = $bam->header;
		$write_alignment = \&write_sam_alignment;
	}
	elsif ($BAM_ADAPTER eq 'hts') {
		$bam = $in->hts_file;
		$header = $bam->header_read;
		$write_alignment = \&write_hts_alignment;
	}
	else {
		die "unrecognized bam adapter $BAM_ADAPTER!";
	}
	
	# output file
	my $out_name = $file;
	$out_name =~ s/\.bam$/.unique.bam/i;
	my $out = Bio::ToolBox::db_helper::write_new_bam_file($out_name) or 
		die "unable to open output bam file $out_name! $!";
		# this uses low level Bio::DB::Bam object
		# using an unexported subroutine imported as necessary depending on bam availability
	$out->header_write($header);
	print " writing $out_name....\n";
	
	# orthologous reads
	my ($out2, $out2_name);
	if ($do_ortho) {
		$out2_name = $file;
		$out2_name =~ s/\.bam$/.orthologous.bam/i;
		my $out = Bio::ToolBox::db_helper::write_new_bam_file($out2_name) or 
			die "unable to open file $out2_name to write! $!\n";
		$out2->header_write($header);
		print " writing $out2_name....\n";
	}
	
	# walk through reads
	my $check = "=$number";
	while (my $a = $bam->read1($header)) {
		my $name = $a->qname;
		next unless exists $alignments{$name};
		if ($alignments{$name} eq $check) {
			&$write_alignment($out, $a, $header);
			delete $alignments{$name};
			$count++;
		}
		elsif ($alignments{$name} eq '=3' and $do_ortho) {
			&$write_alignment($out2, $a, $header);
			$cross_count++;
		}
		elsif (substr($alignments{$name},0,1) ne '=') {
			# an original first file read that never aligned to second file
			# second file reads should never look like this
			&$write_alignment($out, $a, $header);
			delete $alignments{$name};
			$count++;
		}
	}
	
	printf "  %d alignments were written to %s\n", $count, $out_name;
	if ($do_ortho) {
		printf "  %d %s alignments were written to %s\n", $cross_count, 
			$do_pick ? 'equally mapped' : 'orthologous', $out2_name; 
	}
}

sub get_low_as_metric { 
	return -1 * ($_[0]->aux_get('AS') || 0);
}

sub get_high_as_metric { 
	return $_[0]->aux_get('AS') || 0;
}

sub get_mismatch_metric {
	my $nm = $_[0]->aux_get('NM') || $_[0]->aux_get('nM') || 0;
	my $cg = $_[0]->n_cigar || 0;
	return $nm + $cg;
}

sub write_sam_alignment {
	# wrapper for writing Bio::DB::Sam alignments
	# pass bam, alignment, header
	return $_[0]->write1($_[1]);
}

sub write_hts_alignment {
	# wrapper for writing Bio::DB::HTS alignments
	# pass bam, alignment, header
	return $_[0]->write1($_[2], $_[1]);
}


=cut
Notes on various aligners and metrics....

Novoalign 
	mapq good, AS lower is better, NM mismatches
STAR
	mapq fake, AS higher is better, nM mismatches
TopHat and bowtie2
	mapq good, negative AS, higher is better, NM mismatches
BWA
	mapq good, AS higher is better, no NM tag
bowtie1
	don't know, same as bowtie2????



