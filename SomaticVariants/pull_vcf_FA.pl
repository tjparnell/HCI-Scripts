#!/usr/bin/env perl

use strict;
# use Getopt::Long;
# use Data::Dumper;
use Bio::ToolBox::Data;

# usage statement
unless (@ARGV) {
	print <<END;
Script to extract the fraction alternate (FA) values from somatic VCF files.
If your VCF doesn't have the sample FA attribute, run update_somaticVCF_attributes.pl.
Usage: $0 <file1.vcf> <file2.vcf> ...
END
	exit;
}

# iterate through input files
foreach (@ARGV) {
	my $file = shift @ARGV;
	print "reading file $file....\n";
	
	my $Input = Bio::ToolBox::Data->new(file => $file) or 
		die "unable to open $file!\n";
	die "file is not VCF!!!\n" unless $Input->vcf;
	
	my $Output = Bio::ToolBox::Data->new(columns => [qw(Chromosome Position Alternate)]);
	
	# add columns for each sample
	my $number = $Input->number_columns - 1;
	for my $i (9 .. $number) {
		$Output->add_column($Input->name($i));
	}	
	
	# iterate through table
	$Input->iterate( sub {
		my $row = shift;
		my $attributes = $row->vcf_attributes;
		my @fa;
		for my $i (9 .. $number) {
			push @fa, $attributes->{$i}{FA} || 0;
		}
		$Output->add_row( [$row->seq_id, $row->start, $row->value(4), @fa] );
	} );
	
	my $success = $Output->save( $Input->path . $Input->basename . '_FA.txt');
	print "  wrote file $success\n" if $success;
}




