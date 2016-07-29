#!/usr/bin/env perl

use strict;
use Getopt::Long;
use Bio::ToolBox::Data 1.40;

unless (@ARGV) {
	print <<END;

A script to filter somatic VCF files based on the FA (Frequency of Alternate) 
attribute. If you have not done so, run the update_somaticVCF_attributes.pl 
script to add the FA attribute to both tumor and normal samples. Default 
values will pass everything.

Usage: $0 -i <input.vcf> 
Options: 
  --input <input.vcf>   input vcf file, required
  --out <output.vcf>    default is to overwrite input!!!!
  --fail <failed.vcf>   write out failed vcf file, default null
  --norm <NAME>         Normal sample name, default NORMAL
  --tumor <NAME>        Tumor sample name, default TUMOR
  --nmax <float>        Normal FA max, default 1.0
  --nmin <float>        Normal FA min, default 0
  --tmax <float>        Tumor FA max, default 1.0
  --tmin <float>        Tumor FA min, default 0
  --delta <float>       Minimum delta (Tumor - Normal), default 0
END
	exit;
}

my $file;
my $outfile;
my $failfile;
my $normName = 'NORMAL';
my $tumorName = 'TUMOR';
my $normMin  = 0;
my $normMax  = 1.0;
my $tumorMin = 0;
my $tumorMax = 1.0;
my $delta = 0;

GetOptions( 
	'input=s'       => \$file, # input vcf
	'output=s'      => \$outfile, # output vcf
	'fail=s'        => \$failfile, # failed output vcf
	'norm=s'        => \$normName, # normal name
	'tumor=s'       => \$tumorName, # tumor name
	'nmax=f'        => \$normMax, 
	'nmin=f'        => \$normMin,
	'tmax=f'        => \$tumorMax, 
	'tmin=f'        => \$tumorMin,
	'delta=f'       => \$delta,
) or die "bad options!\n";

$outfile ||= $file;
$outfile =~ s/\.gz$//i; # don't mess with bgzip files

# Open VCF file
my $Data = Bio::ToolBox::Data->new(file => $file) or 
	die "unable to open $file!\n";
die "file is not VCF!!!\n" unless $Data->vcf;

# Prepare failed VCF file
my $Fail;
if ($failfile) {
	$failfile =~ s/\.gz$//; # don't mess with bgzip files
	$Fail = $Data->duplicate;
	my $head = $Fail->vcf_headers;
	$head->{FILTER}{FADelta} = 
		qq(ID=FADelta,Description="Difference between Tumor and Normal alternate read frequency (FA) below $delta")
		if $delta;
	$head->{FILTER}{NormFAMin} =
		qq(ID=NormFAMin,Description="Normal alternate read frequency (FA) below $normMin")
		if $normMin;
	$head->{FILTER}{NormFAMax} = 
		qq(ID=NormFAMax,Description="Normal alternate read frequency (FA) exceeds $normMax")
		if $normMax;
	$head->{FILTER}{TumorFAMin} = 
		qq(ID=TumorFAMin,Description="Tumor alternate read frequency (FA) below $tumorMin")
		if $tumorMin;
	$head->{FILTER}{TumorFAMax} = 
		qq(ID=TumorFAMax,Description="Tumor alternate read frequency (FA) exceeds $tumorMax")
		if $tumorMax;
	$Fail->rewrite_vcf_headers;
}

# Filter
my @tosses;
$Data->iterate( \&filter_on_FA );
printf " filtered %s variants, kept %s variants\n", 
	scalar(@tosses), ($Data->last_row - scalar(@tosses));
$Data->delete_row(@tosses);

# Write
my $s = $Data->save($outfile) or die "unable to write $outfile!\n";
print " wrote file $s\n" if $s;
if ($Fail) {
	my $f = $Fail->save($failfile);
	print " wrote failed file $f\n";
}









# filter subroutine
sub filter_on_FA {
	my $row = shift;
	my $att = $row->vcf_attributes;
	return unless exists $att->{$normName}{FA};
	return unless exists $att->{$tumorName}{FA};
	my $check = 0;
	my $filter;
	if ($att->{$normName}{FA} < $normMin) {
		$check++;
		$filter = 'NormFAMin';
	}
	if ($att->{$normName}{FA} > $normMax) {
		$check++;
		$filter = $filter ? "$filter;NormFAMax" : 'NormFAMax';
	}
	if ($att->{$tumorName}{FA} < $tumorMin) {
		$check++;
		$filter = $filter ? "$filter;TumorFAMin" : 'TumorFAMin';
	}
	if ($att->{$tumorName}{FA} > $tumorMax) {
		$check++;
		$filter = $filter ? "$filter;TumorFAMax" : 'TumorFAMax';
	}
	if ( ($att->{$tumorName}{FA} - $att->{$normName}{FA}) < $delta ) {
		$check++;
		$filter = $filter ? "$filter;FADelta" : 'FADelta';
	}
	push @tosses, $row->row_index if $check;
	if ($check and $Fail) {
		if ($row->value(6) eq '.' or $row->value(6) eq 'PASS') {
			$row->value(6, $filter);
		}
		else {
			my $v = $row->value(6);
			$v .= ";$filter";
			$row->value(6, $v);
		}
		$Fail->add_row($row);
	}
}



