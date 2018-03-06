#!/usr/bin/perl

use strict;
use Bio::ToolBox::Data::Stream;

unless (@ARGV) {
	print <<HELP;

Add a variant frequency format field to a VCF file. This is 
calculated from the AD tag as alt / (ref+alt). It is stored 
as Format code VF for "Variant Frequency", to distinguish it 
from the Info field AF "Allele Frequency", which is a population 
frequency.

Usage: $0 <input.vcf> <output.vcf>

HELP
	exit;
}

my $file = shift @ARGV;
my $outfile = shift @ARGV || undef;


# Input file
my $Stream = Bio::ToolBox::Data::Stream->new(
	file      => $file,
) or die "unable to open file '$file'!";
die "not a vcf file!" unless $Stream->vcf;


# Output file
unless ($outfile) {
	$outfile = $Stream->path . $Stream->basename . '.VF' . $Stream->extension;
}
my $Out = $Stream->duplicate($outfile);



# Add format to the metadata
my $md = $Out->vcf_headers;
$md->{FORMAT}{VF} = q[<ID=VF,Number=A,Type=Float,Description="Variant Frequency, alt/(ref+alt)">];
$Out->rewrite_vcf_headers;


# Stream and process
my $number = $Stream->last_column;
while (my $row = $Stream->next_row) {
	my $attrib = $row->vcf_attributes;
	foreach my $i (9 .. $number) {
		# allele depths
		my $ad = $attrib->{$i}{AD} || undef;
		unless ($ad) {
			# store a dummy value
			$attrib->{$i}{VF} = 0;
			next;
		}
		
		# calculate frequency, be mindful of multiple alternate depths
		my ($ref, @alts) = split /,/, $ad;
		# my @freqs = map { $_/($_ + $ref) } @alts;
		my @freqs;
		foreach my $a (@alts) {
			# handle exceptions where both ref and alt is 0
			my $denominator = $ref + $a;
			push @freqs, $denominator ? $a / $denominator : 0;
		}
		
		# put back in
		$attrib->{$i}{VF} = join(',', map { sprintf("%.2f", $_) } @freqs);
	}
	
	# write back vcf changes
	$row->rewrite_vcf_attributes;
	
	# write to output
	$Out->write_row($row);
}

$Out->close_fh;
$Stream->close_fh;
print "wrote $outfile\n";



