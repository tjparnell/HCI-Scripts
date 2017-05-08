#!/usr/bin/env perl

# Timothy J. Parnell, PhD
# Huntsman Cancer Institute
# University of Utah
# Salt Lake City, UT 84112
#  
# This package is free software; you can redistribute it and/or modify
# it under the terms of the Artistic License 2.0.  
# 
# Updated versions of this file may be found in the repository
# https://github.com/tjparnell/HCI-Scripts/

use strict;
use Getopt::Long;
use Bio::ToolBox::Data;

my $outfile;
my $index;

GetOptions( 
	'out=s'       => \$outfile,
	'index=i'     => \$index,
) or die "something's wrong with your options! check your command\n";

unless (@ARGV) {
	die " Need to provide input files!\n";
}

my $OutData = Bio::ToolBox::Data->new( columns => [qw(
	Sample Total nonsynonymous_SNV stoploss stopgain frameshift nonframeshift_indel)]
) or die " failed to create Data object!\n";

# walk through each input file
foreach my $file (@ARGV) {
	my $Data = Bio::ToolBox::Data->new(file => $file) or 
		die " unable to load file '$file'!\n";
	my $name = $Data->name($index);
	my $variant_index = $Data->find_column('Consequence') || 
		$Data->find_column('RefSeqVarType') || $Data->find_column('EnsemblVarType');
	die " unable to find the Consequence column!\n" unless defined $variant_index;
	my $nonsyn = 0;
	my $stoploss = 0;
	my $stopgain = 0;
	my $fs = 0;
	my $nonfs = 0;
	my $other = 0;
	$Data->iterate( sub {
		my $row = shift;
		my $cons = $row->value($variant_index);
		next unless $cons;
		if ($cons =~ /nonsynonymous/i or $cons =~ /missense.?variant/i) {
			$nonsyn++;
		}
		$stoploss++ if ($cons =~ /stop.?loss/i);
		$stopgain++ if ($cons =~ /stop.?gain/i);
		if (
			$cons =~ /non.?frame.?shift/i or 
			$cons =~ /in.?frame.?deletion/i or 
			$cons =~ /in.?frame.?insertion/i
		) {
			$nonfs++ ;
		}
		elsif ($cons =~ /frame.?shift/i) {
			$fs++;
		}
	} );
	
	# add to the output table
	my $total = $nonsyn + $stoploss + $stopgain + $fs + $nonfs;
	$OutData->add_row( [$name, $total, $nonsyn, $stoploss, $stopgain, $fs, $nonfs] );
}

my $s = $OutData->save($outfile);
print "wrote file $s\n";

