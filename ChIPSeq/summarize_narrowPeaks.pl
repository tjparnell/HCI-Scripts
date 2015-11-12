#!/usr/bin/perl

use strict;
use Bio::ToolBox::Data;

print "Summarize scores and lengths from narrowPeak files\n";
unless (scalar @ARGV > 2) {
	print " usage: $0 <outfile> <file1.narrowPeak> <file2.narrowPeak> ....\n";
	exit;
}
my $outfile = shift;

# input files
my @inputs;
foreach my $f (@ARGV) {
	print "  loading $f...\n";
	my $Stream = Bio::ToolBox::Data->new(in => $f, stream => 1);
	next unless $Stream;
	my $D = Bio::ToolBox::Data->new();
	$D->add_column($Stream->basename . '_qvalue');
	$D->add_column($Stream->basename . '_length');
	while (my $row = $Stream->next_row) {
		$D->add_row( [$row->value(8), $row->length] );
	}
	push @inputs, $D;
}

# get the biggest table size
my $max = $inputs[0]->last_row;
foreach (@inputs) {
	$max = $_->last_row if $_->last_row > $max;
}
print " maximum table length is $max\n";

# output data table
my $Data = Bio::ToolBox::Data->new();
foreach my $D (@inputs) {
	printf "  adjusting table from %s to $max\n", $D->last_row;
	while ($D->last_row < $max) {
		$D->add_row;
	}
	my $qvalues = $D->column_values(0);
	$Data->add_column($qvalues);
	my $lengths = $D->column_values(1);
	$Data->add_column($lengths);
}
my $w = $Data->save($outfile);
print " wrote $w\n";


