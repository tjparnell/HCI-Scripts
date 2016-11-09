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
use Bio::DB::Sam;

my $VERSION = 1.2;
# version 1.0 used mapping quality score and number of cigar operations
# version 1.1 uses alignment score and sum of mismatches and cigar operations
# version 1.2 improves text output to make counts less confusing

my $description = <<DESCRIPTION;

Picks the best alignment between two species.

Sequence reads from a mixed species sample, e.g. human xenografts in mouse, 
are aligned to each species index. Unique species alignments are kept, while 
alignments matching to both species are either discarded, written to a cross 
species file, or assigned to one of the two species.

Provide two bam files representing alignments to each species index. The 
bam files should have unique alignments only; only the first occurence of 
multiple alignments will be retained. Bam files do not need to be sorted, 
but any sort order will be maintained. Unmapped reads are ignored and 
discarded.

An option is available to pick the best alignment based on one of two 
criteria: the alignment score (AS tag, lower is better), or the sum of 
the number of mismatches (NM tag) and the number of CIGAR operations 
(including insertions, deletions, soft and hard trims). Alignments 
equal in all scores are considered orthologous.

Othologous reads mapping to both indices can either be discarded (default) 
or written to a second bam file. When the best pick option is enabled, 
the equivalently mapping files are written to the orthologous bam file.

Ouput bam files are appended with the extension .unique.bam and 
.orthologous.bam.

DESCRIPTION
my $usage = <<USAGE;

Usage: $0 [options] file1.bam file2.bam 

Options:
	--ortho     Write a second bam file with all orthologous or 
	            equally mapping alignments
	--pick      Assign the best alignment to one species based on 
	            alignment metrics

USAGE


### Options
	# ugh, GetOptions doesn't like it when no options are given
	# so gotta check if we have options or not
my $do_ortho = 0;
my $do_pick = 0;
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
	) or die "bad options!\n$usage";
	($file1, $file2) = @ARGV;
}
else {
	die "bad options!\n$usage";
	exit 1;
}

### Global values
my %alignments;
my $count1best = 0;
my $count2best = 0;
my $countequal = 0;
my $countortho = 0;


### Reads
my $file1count = read_first_bam();
my $file2count = read_second_bam();

# report numbers
printf " There were %d mapped alignments in $file1\n", $file1count;
printf "   %d alignments mapped better than the second file\n", $count1best 
	if $do_pick;
printf " There were %d mapped alignments in $file2\n", $file2count;
printf "   %d alignments mapped better than the first file\n", $count2best 
	;
printf " There were %d equally mapped alignments\n", $countequal if $do_pick;

if ($do_ortho and not $do_pick) {
	printf " There were %d orthologous alignments\n", $countortho;
}


### Write out the reads
write_bam_files($file1, 1);
write_bam_files($file2, 2);



sub read_first_bam {
	# this reads the first bam file, and stores the alignment score (AS tag) and
	# the sum of CIGAR of operations and number of mismatches
	# (mismatches aren't always stored in the CIGAR) for each aligned read
	
	my $in = Bio::DB::Bam->open($file1) or die " Cannot open input Bam file!\n";
	print " Reading $file1...\n";
	my $count = 0;
	
	# header
	my $h = $in->header;
	
	# walk through reads
	while (my $a = $in->read1) {
		next if $a->unmapped;
		$alignments{$a->qname} = sprintf "%d,%d", $a->aux_get('AS'), 
			$a->aux_get('NM') + $a->n_cigar;
# 		printf("   %s has quality %d CIGAR %s, %d operations, alignment score %d, and %d mismatches\n", $a->qname, $a->qual, $a->cigar_str, $a->n_cigar, $a->aux_get('AS'), $a->aux_get('NM'));
		$count++;
# 		exit if $count > 100;
	}
	return $count;
}

sub read_second_bam {
	# this reads the second bam file, and compares the quality and cigar 
	# number with every mapped alignment that matches from the first file
	# alignments are assigned a single value 
	# 1 is keep first file read
	# 2 is keep second file read
	# 3 is orthologous read in both, keep or discard as requested
	# 0 is discard (multiple hit alignment)
	
	my $in = Bio::DB::Bam->open($file2) or die " Cannot open input Bam file!\n";
	print " Reading $file2...\n";
	my $count = 0;

	# header
	my $h = $in->header;
	
	# walk through reads
	while (my $a = $in->read1) {
		next if $a->unmapped;
		my $name = $a->qname;
		$count++;
	
		# check for existing
		if (exists $alignments{$name}) {
			next if length($alignments{$name}) == 1; # a duplicate alignment!
				# we already processed this so move one
		
			# compare
			my ($score1, $error1) = split ',', $alignments{$name};
			if ($do_pick) {
				my $score = $a->aux_get('AS');
				if ($score < $score1) {
					# second alignment is better
					$alignments{$name} = 2;
					$count2best++;
				}
				elsif ($score > $score1) {
					# first alignment is better
					$alignments{$name} = 1;
					$count1best++;
				}
				elsif ($score == $score1) {
					# equal alignment scores, check number of errors 
					my $error = $a->aux_get('NM') + $a->n_cigar;
					
					if ($error < $error1) {
						# second alignment is better
						$alignments{$name} = 2;
						$count2best++;
					}
					elsif ($error > $error1) {
						# first alignment is better
						$alignments{$name} = 1;
						$count1best++;
					}
					else {
						# alignments are equal
						if ($do_ortho) {
							# keep for third orthologous bam file
							$alignments{$name} = 3;
						}
						else {
							# drop the alignment
							$alignments{$name} = 0;
						}
						$countequal++;
					}
				}
				next;
			}
			if ($do_ortho) {
				# keep for third orthologous bam file
				$alignments{$name} = 3;
				$countortho++;
			}
			else {
				$alignments{$name} = 0 if exists $alignments{$name};
			}
		}
		else {
			# great, unique to the second file
			$alignments{$name} = 2;
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
	my $in = Bio::DB::Bam->open($file);
	
	# output file
	my $out_name = $file;
	$out_name =~ s/\.bam$/.unique.bam/i;
	my $h = $in->header;
	my $out = Bio::DB::Bam->open($out_name, 'w') or 
		die "unable to open file $out_name to write! $!\n";
	$out->header_write($h);
	print " writing $out_name....\n";
	
	# orthologous reads
	my ($out2, $out2_name);
	if ($do_ortho) {
		$out2_name = $file;
		$out2_name =~ s/\.bam$/.orthologous.bam/i;
		my $h = $in->header;
		$out2 = Bio::DB::Bam->open($out2_name, 'w') or 
			die "unable to open file $out2_name to write! $!\n";
		$out2->header_write($h);
		print " writing $out2_name....\n";
	}
	
	# walk through reads
	while (my $a = $in->read1) {
		my $name = $a->qname;
		next unless exists $alignments{$name};
		if ($alignments{$name} == $number) {
			$out->write1($a);
			delete $alignments{$name};
			$count++;
		}
		elsif ($alignments{$name} == 3 and $do_ortho) {
			$out2->write1($a);
			$cross_count++;
		}
		elsif (length($alignments{$name}) > 1) {
			# an original first file read that never aligned to second file
			# second file reads never look like this
			$out->write1($a);
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


