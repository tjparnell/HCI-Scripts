#!/usr/bin/env perl

use strict;
use IO::File;
use File::Basename qw(fileparse);

unless (@ARGV) {
	print <<USAGE;
 
  A script to check and fix aberrant base quality scores in Fastq
  files. It checks for scores outside of the quality range of 0-40, or
  ASCII range 33-73, assuming standard Illumina 1.8+ or Sanger fastq
  offset of 33. It will reset any score above 37 to 37, or F in ASCII.
  It will also reset any score below the acceptable range, just in case.
 
  Usage: $0 file_R1.fastq.gz [file_R2.fastq.gz]
  
  One (single-end) or two (paired-end) files can be provided. Files 
  must be gzipped, and outputs will be gzipped.
  
  Output file names will be the same as the input file name(s) but 
  the basename will be appended with '_fixed'. Summary counts of 
  fixed alignments are reported to standard output.
 
USAGE
	exit(0);
}

### Maximum allowed quality score, ASCII score
# see table below
my $max = 70; # quality score 36 or ASCII F


### Parameters
my ($file1, $file2) = @ARGV;


# Illumina table for getting the quality scores
# https://support.illumina.com/help/BaseSpace_OLH_009008/Content/Source/Informatics/BS/QualityScoreEncoding_swBS.htm


### File handles
my ($infh1, $infh2, $outfh1, $outfh2);

# first file
$infh1 = IO::File->new("gzip -dc $file1 |") or die "unable to read '$file1' $!\n";
$outfh1 = open_output($file1);

# second file
if ($file2) {
	$infh2 = IO::File->new("gzip -dc $file2 |") or die "unable to read '$file2' $!\n";
	$outfh2 = open_output($file2);
}



### Search
my $bad_score = 0;
my $bad_read = 0;
my $total = 0;
while (my $line1 = $infh1->getline) {
	# remaining lines from file 1 for read
	my $line2 = $infh1->getline or die "bad fastq formatting!\n";
	my $line3 = $infh1->getline or die "bad fastq formatting!\n";
	my $line4 = $infh1->getline or die "bad fastq formatting!\n";
	
	# process
	$line4 = process_quality($line4);
	$outfh1->print($line1 . $line2 . $line3 . $line4);
	
	# second file
	if ($file2) {
		my $line21 = $infh2->getline or die "bad fastq formatting!\n";
		my $line22 = $infh2->getline or die "bad fastq formatting!\n";
		my $line23 = $infh2->getline or die "bad fastq formatting!\n";
		my $line24 = $infh2->getline or die "bad fastq formatting!\n";
		
		# process
		$line24 = process_quality($line24);
		$outfh2->print($line21 . $line22 . $line23 . $line24);
	}
	
	$total++;
}

### Done
$infh1->close;
$infh2->close if $infh2;
$outfh1->close;
$outfh2->close if $outfh2;
print " Processed $total reads\n";
print "  There were $bad_score out-of-bound quality scores in $bad_read reads\n";







sub open_output {
	my $file = shift;
	my ($basename, $path, $extension) = fileparse($file, qw(.txt.gz .fastq.gz .fq.gz));
	my $filename = $path . $basename . '_fixed' . $extension;
	print " Writing to $filename\n";
	my $fh = IO::File->new("| gzip >$filename") or 
		die "cannot write to compressed file '$filename' $!\n";
	return $fh;
}


sub process_quality {
	my $line = shift;
	chomp $line;
	
	# process scores
	my @scores = map { ord($_) } split(q(), $line);
	my @new;
	my $error = 0;
	foreach my $s (@scores) {
		if ($s > $max) {
			# out-of-bounds value
			$error++;
			$s = $max; # set to maximum allowed
		}
		elsif ($s < 33) {
			# ASCII control character
			$error++;
			$s = 33; # sets to ! in ASCII
		}
		push @new, $s;
	}
	
	# update counts
	$bad_score += $error;
	$bad_read++ if $error;
	return sprintf("%s\n", join(q(), map { chr($_) } @new));
}

=cut

Illumina table for getting the quality scores
https://support.illumina.com/help/BaseSpace_OLH_009008/Content/Source/Informatics/BS/QualityScoreEncoding_swBS.htm

Character	ASCII	Quality
!	33	0
"	34	1
#	35	2
$	36	3
%	37	4
&	38	5
'	39	6
(	40	7
)	41	8
*	42	9
+	43	10
,	44	11
-	45	12
.	46	13
/	47	14
0	48	15
1	49	16
2	50	17
3	51	18
4	52	19
5	53	20
6	54	21
7	55	22
8	56	23
9	57	24
:	58	25
;	59	26
<	60	27
=	61	28
>	62	29
?	63	30
@	64	31
A	65	32
B	66	33
C	67	34
D	68	35
E	69	36
F	70	37
G	71	38
H	72	39
I	73	40

