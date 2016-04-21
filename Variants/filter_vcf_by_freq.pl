#!/usr/bin/perl

use strict;
use Bio::ToolBox::Data::Stream;

unless (@ARGV) {
	print <<HELP;

filter a single sample VCF file by it's frequency of alt/ref

Usage: $0 <FREQ> <input.vcf> <output.vcf>

FREQ should be a floating point value between 0 and 1
ouput.vcf is optional, default is input.filter.vcf

HELP
	exit;
}

my $freq = shift @ARGV;
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
	my $ad = $attrib->{9}{AD} || undef;
	if ($ad) {
		my ($ref, $alt) = split /,/, $ad;
		if ( ($alt/$ref) >= $freq) {
			$Out->write_row($row);
			$passCount++;
		}
		else {
			$failCount++;
		}
	}
	else {
		$failCount++;
	}
}

print "passed $passCount variants, failed $failCount variants\n";
