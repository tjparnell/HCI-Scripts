#!/usr/bin/perl

use strict;
use Bio::ToolBox::Data::Stream;

unless (@ARGV) {
	print <<HELP;

Filter a single sample VCF file by it's frequency of alternate reads. 
Keeps all variants that are equal or greater to the given frequency threshold.

Usage: $0 <FREQ> <input.vcf> <output.vcf>

FREQ should be a floating point value between 0 and 1
ouput.vcf is optional, default is input.filter.vcf

HELP
	exit;
}

my $cutoff = shift @ARGV;
my $file = shift @ARGV;
my $outfile = shift @ARGV;

my $Stream = Bio::ToolBox::Data::Stream->new(
	file      => $file,
) or die "unable to open file '$file'!";

die "not a vcf file!" unless $Stream->vcf;

unless ($outfile) {
	$outfile = $Stream->path . $Stream->basename . '.filter' . $Stream->extension;
}


my $Out = $Stream->duplicate($outfile);

my $passCount = 0;
my $failCount = 0;
while (my $row = $Stream->next_row) {
	my $attrib = $row->vcf_attributes;
	my $freq = $attrib->{9}{FREQ} || $attrib->{9}{FA} || undef;
	if ($freq =~ /%$/) {
		# reported as actual percentage, convert to floating point value
		$freq =~ s/%$//;
		$freq /= 100;
	}
	if (not defined $freq) {
		my $ad = $attrib->{9}{AD} || undef;
		my ($ref, $alt) = split /,/, $ad;
		$freq = $ref ? ($alt/($alt + $ref)) : 0;
	}
	if ( $freq >= $cutoff) {
		$Out->write_row($row);
		$passCount++;
	}
	else {
		$failCount++;
	}
}

print "passed $passCount variants, failed $failCount variants\n";
