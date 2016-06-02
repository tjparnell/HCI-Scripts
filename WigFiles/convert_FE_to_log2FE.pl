#!/usr/bin/perl

use strict;
use IO::File;
use Bio::ToolBox::big_helper qw(wig_to_bigwig_conversion);
use constant LOG2 => log(2);
use constant LOG10 => log(10);

print "\n Script to convert Macs2 FE files to log2 FE files\n";

# files
my $infile = shift @ARGV or 
	die "usage: $0 file_FE\n";


# using logFE or FE
# assuming the input file is appropriately named
my $delog10;
if ($infile =~ /logFE\.\w{2,}$/) {
	$delog10 = 1;
	print "converting from log10 FE\n";
}
else {
	print "did not detect logFE in file name. Assuming linear FE\n";
}

# convert to bedgraph as necessary
my $bdgfile;
my $use_bw = 0; # flag to indicate we are working with bigwig files

if ($infile =~ /\.bw$/) {
	$use_bw = 1;
	$bdgfile = $infile;
	$bdgfile =~ s/bw$/bdg/;
	system('bigWigToBedGraph', $infile, $bdgfile) == 0 or 
		die "unable to execute bigWigToBedGraph!\n";
}
else {
	$bdgfile = $infile;
}


# open file handles
my $newfile = $bdgfile;
$newfile =~ s/(log)?FE\.bdg$/log2FE.bdg/;
my $infh = IO::File->new($bdgfile, 'r') or 
	die "unable to open input bedGraph file $bdgfile!\n";
my $outfh = IO::File->new($newfile, 'w') or 
	die "unable to open output bedGraph file $newfile!\n";


# convert scores
if ($delog10) {
	# from logFE to log2
	while (my $line = $infh->getline) {
		chomp $line;
		my @data = split /\t/, $line;
		if ($data[3] =~ /^\-?\d+\.?\d*$/) {
			$data[3] = log(10 ** $data[3]) / LOG2;
		}
		$outfh->print( join("\t", @data) . "\n");
	}
}
else {
	# convert scores from linear FE to log2
	while (my $line = $infh->getline) {
		chomp $line;
		my @data = split /\t/, $line;
		if ($data[3] =~ /^\-?\d+\.?\d*$/ and $data[3] != 0) {
			$data[3] = log($data[3]) / LOG2;
		}
		$outfh->print( join("\t", @data) . "\n");
	}
}


# finish
$infh->close;
$outfh->close;

# convert back to bigwig
if ($use_bw) {
	my $newbw = wig_to_bigwig_conversion(
		wig     => $newfile,
		db      => $infile,
	) or die "failed to convert bedgraph to bigwig!\n";
	unlink($bdgfile, $newfile);
	print "converted $infile to $newbw\n";
}
else {
	print "converted $bdgfile to $newfile\n";
}



