#!/usr/bin/env perl

# pull reads from fastq
# this requires the library from UMIScripts

use strict;
use IO::File;
use List::Util qw(sum0);
use Bio::UMIScripts::FastqHelper;

unless (@ARGV) {
	print "\nA script to pull out Fastq reads by name\n\n$0 <out> <list.txt> <fastq1> ...\n\n";
	exit;
}

my $outfile = shift @ARGV;
my $listfile = shift @ARGV;

# get list
my $fh = IO::File->new($listfile) or die "unable to open $listfile! $!\n";
my %list;
while (my $line = $fh->getline) {
	chomp $line;
	$list{$line} = 0;
}
$fh->close;
printf " collected %d names from list\n", scalar(keys %list);

# open output
my $outfh = write_filehandle($outfile);


# read files
while (@ARGV) {
	my $infile = shift @ARGV;
	my $infh = read_fastq_filehandle($infile);
	print " reading $infile...\n";
	while (my $read = get_fastq_read($infh)) {
		my $n = substr($read->[NAME],1);
		if (exists $list{$n}) {
			$outfh->print( fastq_string($read) );
			$list{$n}++;
		}
	}
	$infh->close;
}
$outfh->close;


# results
printf " pulled %d reads\n", sum0(values %list);
my @nocounts = grep {$list{$_} == 0} keys %list;
if (@nocounts) {
	printf " %d reads were not found\n", scalar(@nocounts);
}


