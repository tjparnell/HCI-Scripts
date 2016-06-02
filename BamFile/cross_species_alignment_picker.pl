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


my $description = <<DESCRIPTION;

Picks the best alignment between two species.

Sequence reads from a mixed species sample, e.g. human xenografts in mouse, 
are aligned to each species index. Unique species alignments are kept, while 
alignments matching to both species are either discarded, written to a cross 
species file, or assigned to one of the two species.

Provide two bam files representing alignments to each species index. The 
bam files should have unique alignments only; multiple alignments will 
probably be ignored. Bam files do not need to be sorted, but any sort order 
will be maintained. Unmapped reads are ignored and discarded.

An option is available to pick the best alignment based on the mapping 
quality score and the number of CIGAR operations (on the assumption that 
a simpler alignment strategy with fewest number of mismatches and indels
is better). 

Othologous reads mapping to both indices can either be discarded (default) 
or written to a second bam file. When the best pick option is enabled, 
the equivalently mapping files are written to the orthologous bam file.

Ouput bam files are appended with the extension .unique.bam and 
.orthologous.bam.

DESCRIPTION
my $usage = <<USAGE;

Usage: $0 [options] file1.bam file2.bam 

Options:
	--ortho     Write a second bam file with all orthologous reads
	--pick      Assign the best alignment to one species based on 
	            mapping score and the fewest number of CIGAR operations
USAGE


### Options
	# ugh, GetOptions doesn't like it when no options are given
	# so gotta check if we have options or not
my $do_ortho = 0;
my $do_pick = 0;
my ($file1, $file2);
if (scalar @ARGV == 0) {
	print $description $usage;
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
read_first_bam();
read_second_bam();

# report numbers
if ($do_pick) {
	printf " There were %d alignments best aligned in file 1\n", $count1best;
	printf " There were %d alignments best aligned in file 2\n", $count2best;
	printf " There were %d equally aligned alignments\n", $countequal;
}
if ($do_ortho) {
	printf " There were %d orthologous alignments\n", $countortho;
}


### Write out the reads
write_bam_files($file1, 1);
write_bam_files($file2, 2);



sub read_first_bam {
	# this reads the first bam file, and stores the quality and cigar number
	# for each aligned read
	
	my $in = Bio::DB::Bam->open($file1) or die " Cannot open input Bam file!\n";
	print " Reading $file1...\n";
	my $count = 0;
	
	# header
	my $h = $in->header;
	
	# walk through reads
	while (my $a = $in->read1) {
		next if $a->unmapped;
		$alignments{$a->qname} = sprintf "%d\t%d", $a->qual, $a->n_cigar;
# 		printf("   %s has quality %d and %d CIGAR operations\n", $a->qname, $a->qual, $a->n_cigar) if $count < 10;
		$count++;
	}
	print "  read $count mapped alignments.\n";
}

sub read_second_bam {
	# this reads the second bam file, and compares the quality and cigar 
	# number with every mapped alignment that matches from the first file
	# alignments are assigned a single value 
	# 1 is keep first file read
	# 2 is keep second file read
	# 3 is orthologous read in both, keep or discard as requested
	
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
			my ($qual, $nop) = split '\t', $alignments{$name};
			if ($do_pick) {
				if ($a->qual > $qual) {
					# second alignment is better
					$alignments{$name} = 2;
					$count2best++;
				}
				elsif ($a->qual < $qual) {
					# first alignment is better
					$alignments{$name} = 1;
					$count1best++;
				}
				elsif ( ($a->qual == $qual) and ($a->n_cigar < $nop) ) {
					# second alignment is better
					$alignments{$name} = 2;
					$count2best++;
				}
				elsif ( ($a->qual == $qual) and ($a->n_cigar > $nop) ) {
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
						delete $alignments{$name};
					}
					$countequal++;
				}
				next;
			}
			if ($do_ortho) {
				# keep for third orthologous bam file
				$alignments{$name} = 3;
				$countortho++;
			}
			else {
				delete $alignments{$name} if exists $alignments{$name};
			}
		}
		else {
			# great, unique to the second file
			$alignments{$name} = 2;
		}
	}
	print "  read $count mapped alignments.\n";
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
	printf("  %d orthologous alignments were written to %s\n", $cross_count, 
		$out2_name) if $do_ortho; 
}


