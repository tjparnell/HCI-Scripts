#!/usr/bin/perl -w

use strict;
use IO::File;
use File::Basename qw(fileparse);

unless (@ARGV) {
	print <<USAGE;
 
  A script to pull out the remaining un-aligned fastq reads.
  Find the last alignment read name by running 
    "samtools view file.raw.bam | tail"
  The raw bam is very likely truncated and samtools will give an error message.
  Record the name of of the last alignment
 
  Usage: $0 'name' file1.txt.gz file2.txt.gz
  
  Where 'name' is the read name of the last alignment. 
  Output file names will be the same as the input file name but the basename 
  will be appended with '_remaining'.
 
USAGE
	exit(0);
}

### Parameters
my ($name, $file1, $file2) = @ARGV;




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
my $found = 0;
my $align_count = 1;
my $todo_count = 1;

while (my $line1 = $infh1->getline) {
	# remaining lines from file 1
	my $line2 = $infh1->getline;
	my $line3 = $infh1->getline;
	my $line4 = $infh1->getline;
	
	# second file
	my ($line21, $line22, $line23, $line24);
	if ($file2) {
		$line21 = $infh2->getline;
		$line22 = $infh2->getline;
		$line23 = $infh2->getline;
		$line24 = $infh2->getline;
	}
	
	# write if we found what we're looking for
	if ($found) {
		$outfh1->print($line1 . $line2 . $line3 . $line4);
		if ($file2) {
			$outfh2->print($line21 . $line22 . $line23 . $line24);
		}
		$todo_count++;
	}
	elsif ($line1 =~ /$name/) {
		# we got it!!!!!
		# start writing the next round
		$found = 1;
		$align_count++;
		print "found match at read $align_count\n";
	}
	else {
		# drat! nothing yet
		$align_count++;
	}
}

### Done
$infh1->close;
$infh2->close if $infh2;
$outfh1->close;
$outfh2->close if $outfh2;
print " Skipped $align_count reads\n";
print " Wrote $todo_count reads\n";



sub open_output {
	my $file = shift;
	my ($basename, $path, $extension) = fileparse($file, qw(.txt.gz .fastq.gz .fq.gz));
	my $filename = $path . $basename . '_remaining' . $extension;
	if ($basename =~ /^(.+)([_\.\\\/][12])$/) {
		# retain the part number
		$filename = $path. $1 . '_remaining' . $2 . $extension;
	}
	print " Writing to $filename\n";
	my $fh = IO::File->new("| gzip >$filename") or 
		die "cannot write to compressed file '$filename' $!\n";
	return $fh;
}




