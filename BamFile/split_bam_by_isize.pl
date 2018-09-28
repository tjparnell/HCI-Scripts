#!/usr/bin/env perl

# documentation at end of file

# Updated versions of this file may be found in the repository
# https://github.com/tjparnell/HCI-Scripts/

use strict;
use Getopt::Long;
use Pod::Usage;
use Bio::ToolBox::db_helper 1.60 qw(
	open_db_connection 
	low_level_bam_fetch
	$BAM_ADAPTER
);
use Bio::ToolBox::utility;
my $VERSION = '1.62';


print "\n A script to split a paired-end bam file by insert sizes\n\n";

### Quick help
unless (@ARGV) { 
	# when no command line options are present
	# print SYNOPSIS
	pod2usage( {
		'-verbose' => 0, 
		'-exitval' => 1,
	} );
}


### Command line options
my (
	$infile, 
	$outfile, 
	$failfile,
	$minsize, 
	$maxsize, 
	$AT_ends,
	$min_qual,
	$quick,
	$help,
	$print_version,
);
my @size_list;
GetOptions( 
	'in=s'       => \$infile, # the input bam file path
	'out=s'      => \$outfile, # name of output file 
	'fail=s'     => \$failfile, # the name of the failed reads file
	'min=i'      => \$minsize, # the minimum cutoff size for paired-read segments
	'max=i'      => \$maxsize, # the maximum cutoff size for paired-read segments
	'size=s'     => \@size_list, # a list of sizes to select
	'at'         => \$AT_ends, # discard non-AT ends
	'qual=i'     => \$min_qual, # minimum quality
	'quick!'     => \$quick, # don't look for paired reads, go on flags only
	'help'       => \$help, # request help
	'version'    => \$print_version, # print the version
) or die " unrecognized option(s)!! please refer to the help documentation\n\n";

# Print help
if ($help) {
	# print entire POD
	pod2usage( {
		'-verbose' => 2,
		'-exitval' => 1,
	} );
}

# Print version
if ($print_version) {
	print " Biotoolbox script split_bam_by_isize.pl, version $VERSION\n\n";
	exit;
}




### Check for required values and set defaults
# input file
unless ($infile) {
	if (@ARGV) {
		$infile = shift @ARGV;
	}
	else {
		die "  An input BAM file must be specified!\n";
	}
}

# sizes to select
my @sizes; # an array of arrays of the sizes to include
if (@size_list) {
	foreach (@size_list) {
		my @s = split /-/;
		die " Improperly formatted size range [$_]!\n" if @s != 2;
		push @sizes, [ @s ];
	}
}
else {
	unless (defined $minsize) {
		$minsize = 100; # set default to 100 bp
		print " Using default minimum size of 100 bp\n";
	}
	unless (defined $maxsize) {
		$maxsize = 200; # set default to 200 bp
		print " Using default maximum size of 200 bp\n";
	} 
	push @sizes, [ ($minsize, $maxsize) ];
}

# identify the lowest and highest sizes
my $lowest = $sizes[0][0];
my $highest = $sizes[0][1];
foreach (@sizes) {
	$lowest = $_->[0] if $_->[0] < $lowest;
	$highest = $_->[1] if $_->[1] > $highest;
}

unless ($outfile) {
	$outfile = $infile;
	$outfile =~ s/\.bam$//;
	$outfile =~ s/\.sorted//;
}

# minimum quality
$min_qual ||= 0;
	

### Initialize counts
my $read_count           = 0;
my $pair_count           = 0;
my $quality_count        = 0;
my $missing_paired_count = 0; 
my $no_left_mate_count   = 0;
my $non_paired_count     = 0;
my $improper_count       = 0;
my $mate_unmapped_count  = 0;
my $same_strand_count    = 0;
my $diffchromo_count     = 0;
my $non_AT_end_count     = 0;
my $toosmall_count       = 0;
my $toobig_count         = 0;
my $just_right_count     = 0;

# buffer for the left hand paired reads until we come to the right hand read
my %buffer;




### Open BAM files
print " Opening bam files....\n";

# input file
my $in_sam = open_db_connection($infile) 
	or die " unable to open input bam file '$infile'!\n";
print "   input file '$infile'\n";

# input header
my ($header, $write_alignment);
if ($BAM_ADAPTER eq 'sam') {
	$header = $in_sam->bam->header;
	$write_alignment = \&write_sam_alignment;
}
elsif ($BAM_ADAPTER eq 'hts') {
	$header = $in_sam->hts_file->header_read;
	$write_alignment = \&write_hts_alignment;
}

# output files
foreach (@sizes) {
	# we will open a separate file for each size range
	
	# generate specific file name
	my $bam_file = $outfile . '.' . $_->[0] . '_' . $_->[1] . '.bam';
	
	# open bam file
	my $bam = Bio::ToolBox::db_helper::write_new_bam_file($bam_file) or 
		die "unable to open output bam file $bam_file! $!";
		# using an unexported subroutine imported as necessary depending on bam availability
	$bam->header_write($header);
	print "   output file '$bam_file'\n";
	
	# store
	$_->[2] = 0; # a count for the number of pairs written to this file
	$_->[3] = $bam_file;
	$_->[4] = $bam;
}

# failed file
my $failed; # for storing the failed bam file object
if ($failfile){
	$failfile .= '.bam' unless $failfile =~ /\.bam$/;
	$failed = Bio::ToolBox::db_helper::write_new_bam_file($failfile) or 
		die "unable to open failed output bam file '$failfile' for writing!\n";
	$failed->header_write($header);
	print "   output file '$failfile'\n";
}




### Start conversion
print " Splitting reads...\n";
for my $tid (0 .. $in_sam->n_targets - 1) {
	printf "  sequence %s....\n", $in_sam->target_name($tid);
	
	if ($quick) {
		low_level_bam_fetch($in_sam, $tid, 0, $in_sam->target_len($tid), \&quick_callback, 1);
	}
	else {
		low_level_bam_fetch($in_sam, $tid, 0, $in_sam->target_len($tid), \&paired_callback, 1);
	}
	
	# check for orphans
	if (%buffer) {
		$missing_paired_count += scalar(keys %buffer);
		undef %buffer;
	}
}





### Finish up 
print " Finished splitting\n Bam files need to be sorted and indexed\n\n";


# Print summaries
print "\n There were ", format_with_commas($read_count)," total mapped reads\n";
print "\n There were ", format_with_commas($quality_count)," reads that failed" . 
	" minimum quality of $min_qual\n";
print " There were ", format_with_commas($non_paired_count), " non-paired reads\n" if 
	$non_paired_count;
print " There were ", format_with_commas($pair_count), " total alignment pairs\n";

print "   " . format_with_commas($improper_count) . " (". percent_pc($improper_count) . 
	") pairs were improper\n" if $improper_count;
print "     " . format_with_commas($mate_unmapped_count) . " (". percent_pc($mate_unmapped_count) . 
	") pairs had an unmapped mate\n" if $mate_unmapped_count;
print "     " . format_with_commas($same_strand_count) . " (". percent_pc($same_strand_count) . 
	") pairs had mates on the same strand\n" if $same_strand_count;
print "     " . format_with_commas($diffchromo_count) . " (". percent_pc($diffchromo_count) . 
	") pairs had mates on different chromosomes\n" if $diffchromo_count;

print "   " . format_with_commas($toosmall_count + $toobig_count + $just_right_count) . 
	" (" . percent_pc($toosmall_count + $toobig_count + $just_right_count) . 
	") pairs were proper\n";
print "     " . format_with_commas($toosmall_count) . " (". percent_pc($toosmall_count) . 
	") pairs had insertions below the minimum $lowest bp\n";
print "     " . format_with_commas($toobig_count) . " (". percent_pc($toobig_count) . 
	") pairs had insertions above the maxiumum $highest bp\n";
print "     " . format_with_commas($just_right_count) . " (". percent_pc($just_right_count) . 
	") pairs had insertions of acceptable size\n";

my $failed_sum = $non_AT_end_count + $missing_paired_count + $no_left_mate_count;
print "   " . format_with_commas($failed_sum) . " (" . percent_pc($failed_sum) . 
	") pairs failed to write\n" if $failed_sum;
print "     " . format_with_commas($non_AT_end_count) . " (". percent_pc($non_AT_end_count) . 
	") pairs had one or more non-AT ends\n" if $non_AT_end_count;
print "     " . format_with_commas($missing_paired_count) . " (". percent_pc($missing_paired_count) . 
	") pairs had a missing right mate\n" if $missing_paired_count;
print "     " . format_with_commas($no_left_mate_count) . " (". percent_pc($no_left_mate_count) . 
	") pairs had a missing left mate\n" if $no_left_mate_count;


foreach (@sizes) {
	print " " . format_with_commas($_->[2]) . " (" . percent_pc($_->[2]) . 
		") pairs were written to file '$_->[3]'\n";
}
print "\n";






### Subroutines

sub percent_pc {
	# for calculating the percent of pair_count (total)
	my $count = shift;
	return sprintf "%.2f%%", ($count / $pair_count) * 100;
}

sub paired_callback {
	my $a = $_[0];
	$read_count++;
	
	# mapping quality
	if ($a->qual < $min_qual) {
		&$write_alignment($failed, $a) if $failed;
		$quality_count++;
		return;
	}
	
	# check
	unless ($a->paired) {
		$non_paired_count++;
		&$write_alignment($failed, $a) if $failed;
		return;
	}
	unless ($a->proper_pair) {
		
		# determine why
		if ($a->tid != $a->mtid) {
			if ($a->tid < $a->mtid) {
				$pair_count++;
				$improper_count++;
				$diffchromo_count++;
			}
		}
		elsif ($a->munmapped) {
			$pair_count++;
			$improper_count++;
			$mate_unmapped_count++;
		}
		elsif ($a->reversed == $a->mreversed) {
			# this implies both mates are on the same strand 
			if ($a->pos < $a->mpos) {
				$pair_count++;
				$improper_count++;
				$same_strand_count++;
			}
		}
		&$write_alignment($failed, $a) if $failed;
		return;
	}
	
	# Process a right hand mate
	if ($a->reversed) {
		
		# look for its pair and write
		if (exists $buffer{ $a->qname }) {
			# awesome, we now have both mates to write
			
			# check for AT ends
			if ($AT_ends) {
				my $dna1 = $buffer{ $a->qname }->qseq || 'N';
				my $dna2 = $a->qseq || 'N';
				unless ($dna1 =~ /^[AaTt]/ and $dna2 =~ /^[AaTt]/) {
					$non_AT_end_count++;
					&$write_alignment($failed, $a) if $failed;
					&$write_alignment($failed, $buffer{$a->qname}) if $failed;
					delete $buffer{ $a->qname };
					return;
				}
			}
			
			# write the alignments
			write_out_alignments( $buffer{ $a->qname }, $a);
			delete $buffer{ $a->qname };
		}
		else {
			# the left mate should exist, we shouldn't reach this far for 
			# it not be present
			# most likely explanation is the size was not right
			if ($a->isize >= $lowest and $a->isize <= $highest) {
				# the right size but no left mate, that is odd
				$no_left_mate_count++;
				&$write_alignment($failed, $a) if $failed;
			}
		}
		return;
	}
	
	
	# the rest of this sub assume a left hand mate of a pair
	$pair_count++;
	
	# check size
	my $size = $a->isize;
	if ($size < $lowest) {
		$toosmall_count++;
		&$write_alignment($failed, $a) if $failed;
		return;
	}
	elsif ($size > $highest) {
		$toobig_count++;
		&$write_alignment($failed, $a) if $failed;
		return;
	}
	else {
		$just_right_count++;
		
		# store it in the buffer to write once we have the pair
		$buffer{ $a->qname } = $a;
	}
}

sub quick_callback {
	my $a = $_[0];
	$read_count++;
	
	# mapping quality
	if ($a->qual < $min_qual) {
		&$write_alignment($failed, $a) if $failed;
		$quality_count++;
		return;
	}
	
	# check
	unless ($a->paired) {
		$non_paired_count++;
		&$write_alignment($failed, $a) if $failed;
		return;
	}
	unless ($a->proper_pair) {
		$improper_count++;
		
		# determine why
		if ($a->tid != $a->mtid) {
			if ($a->tid < $a->mtid) {
				$pair_count++;
				$improper_count++;
				$diffchromo_count++;
			}
		}
		elsif ($a->munmapped) {
			$pair_count++;
			$improper_count++;
			$mate_unmapped_count++;
		}
		elsif ($a->reversed == $a->mreversed) {
			# this implies both mates are on the same strand 
			if ($a->pos < $a->mpos) {
				$pair_count++;
				$improper_count++;
				$same_strand_count++;
			}
		}
		&$write_alignment($failed, $a) if $failed;
		return;
	}
	$pair_count++ unless ($a->reversed); # count only left reads as pair
	
	# check size
	my $size = $a->isize;
	if ($size < $lowest) {
		$toosmall_count++;
		&$write_alignment($failed, $a) if $failed;
		return;
	}
	elsif ($size > $highest) {
		$toobig_count++;
		&$write_alignment($failed, $a) if $failed;
		return;
	}
	else {
		$just_right_count++;
		# write it immediately
			
		# check for AT ends if requested
		if ($AT_ends) {
			my $dna = $a->qseq || 'N';
			unless ($dna =~ /^[AaTt]/) {
				$non_AT_end_count++;
				&$write_alignment($failed, $a) if $failed;
				return;
			}
		}
		
		# write
		write_out_alignments($a);
	}
}


sub write_out_alignments {
	my ($a1, $a2) = @_;
	my $length = $a1->isize;
	
	# find the right size range and write to the appropriate file
	foreach my $size (@sizes) {
		# check sizes
		if (
			$length >= $size->[0] and
			$length <= $size->[1]
		) {
			&$write_alignment($size->[4], $a1);
			&$write_alignment($size->[4], $a2) if defined $a2;
			$size->[2] += 1 unless ($a1->reversed); # count the left ones only
		}
	}
}

sub write_sam_alignment {
	# wrapper for writing Bio::DB::Sam alignments
	return $_[0]->write1($_[1]);
}

sub write_hts_alignment {
	# wrapper for writing Bio::DB::HTS alignments
	return $_[0]->write1($header, $_[1]);
}




__END__

=head1 NAME

split_bam_by_isize.pl

A script to selectively write out paired-end alignments of a given insertion size.

=head1 SYNOPSIS

split_bam_by_isize.pl [--options] <file.bam>
  
  Options:
  --in <file.bam>
  --min <integer>       (100)
  --max <integer>       (200)
  --size <min-max>
  --out <filename>
  --fail <filename>
  --qual <integer>      0
  --at
  --quick
  --version
  --help

=head1 OPTIONS

The command line flags and descriptions:

=over 4

=item --in <file.bam>

Specify the file name of a binary Bam file of paired-end alignments 
as described for Samtools. Use samtools to convert text Sam files 
to binary Bam files. The file should be sorted and indexed; this 
program should be able to do it for you automatically if it is not 
(assuming the bam directory is writeable).

=item --min <integer>

Optionally specify the minimum size of fragment to include in the output
Bam file. Not required if the --size option is used. The default value 
is 100 bp.

=item --max <integer>

Optionally specify the maximum size of fragment to include in the output
Bam file. Not required if the --size option is used. The default value 
is 200 bp.

=item --size <min-max>

When multiple size ranges are desired, they may be specified using the 
size option. Define the minimum and maximum size as a range separated by 
a dash (no spaces). Use this option repeatedly for multiple size ranges. 
Size ranges may overlap. The size option takes precedence over the --min 
and --max options.

=item --out <filename>

Optionally specify the base name of the output file. The default is to use 
base input name. The output file names are appended with '.$min_$max'. 

=item --fail <filename>

Optionally provide an output file name where failed alignments are stored. 
Here, failed alignments are those that do not pass the insert size 
requirements or are un- or improperly paired. Un-aligned reads are not written. 

=item --qual <integer>

Optionally provide a minimum alignment mapping quality score to accept. 
The default is 0 (accept everything). This filter is applied prior to 
checking pairs, so failure counts may be under estimated due to at least 
one of the alignments failing the mapping quality score.

=item --at

Boolean option to indicate that only alignments whose query sequence 
begins with an [AT] nucleotide should be included in the output Bam file(s). 
Micrococcal nuclease (MNase) cuts (almost?) exclusively at [AT][AT] 
dinucleotides; this option ensures that the fragment is more likely derived 
from a MNase cut.

By default, both mates must pass the AT test prior to writing.

=item --quick

Boolean option to quickly split the Bam file without checking both mates. 
By default, both mates in a pair must exist and be checked prior to writing, 
which increases processing time and memory requirements. Each mate in the 
pair must have the same ID tag to be considered paired.

When both the --quick and --at option (above) are enabled, each read mate is 
checked independently, potentially leading to the possibility of only one 
mate in a pair being written to the output file.

=item --version

Print the version number.

=item --help

Display the POD documentation

=back

=head1 DESCRIPTION

This program will read a BAM file containing paired-read alignments and
write a new BAM file containing only those successfully aligned paired 
reads whose insert size fall within the set minimum and maximum 
lengths inclusively (defaults of 100 and 200 bp, respectively). 

Multiple size ranges may be specified, and pairs within each range are 
written to separate files. 

Additionally, fragments may be checked for the presence of an A/T nucleotide 
at the 5' end of the sequence, as is usually produced by digestion with 
Micrococcal Nuclease. This does NOT ensure that the original cut site was 
AA, TT, or AT, as we only have one half of the site in the read sequence, 
only that it increases the liklihood of being an A/T dinucleotide 
(absolute checking would require looking up the original sequence).

The input and output files are BAM files as described by the Samtools 
project (http://samtools.sourceforge.net). 

The input file should be sorted and indexed prior to sizing. The output 
file(s) will be automatically re-sorted and re-indexed for you.

A number of statistics about the read pairs are also written to standard 
output for your information, including the number of proper alignments 
in each requested size range, the number of alignments that fail the AT 
check (if requested), and the number of improper paired alignments, 
including those whose mates align to the same strand or different 
chromosomes, or pairs with a non-aligned mate. Given the vagaries of 
different alignment possibilities and enumerations, not all possible 
combinations may be accounted. 

If no pairs are written to the output file(s), try enabling the --quick 
option. It may be that only one of the two mates are actually present 
in the Bam file, or that they do not have matching ID numbers.

=head1 AUTHOR

 Timothy J. Parnell, PhD
 Howard Hughes Medical Institute
 Dept of Oncological Sciences
 Huntsman Cancer Institute
 University of Utah
 Salt Lake City, UT, 84112

This package is free software; you can redistribute it and/or modify
it under the terms of the Artistic License 2.0.  
