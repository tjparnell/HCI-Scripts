#!/usr/bin/env perl

use strict;
use Statistics::Lite qw(min max sum mean median stddev);
use Bio::ToolBox::Data;

unless (scalar @ARGV > 1) {
	print "summarize fraction of coverage depth from 1 or more bedGraph files\n";
	print "usage: $0 <outfile> <file1.bdg.gz> <file2.bdg> ...\n";
	exit;
}


my $outfile = shift (@ARGV);
my $Data = Bio::ToolBox::Data->new();
my $depth_i = $Data->add_column('Depth');
for (my $d = 0; $d <= 500; $d++) {
	$Data->add_row([$d]);
}

while (@ARGV) {
	my $file = shift @ARGV;
	print "loading file $file\n";
	my $Stream = Bio::ToolBox::Data->new(
		stream  => 1,
		in      => $file,
	);
	unless ($Stream) {
		warn "bad $file $!";
		next;
	}
	
	# get counts
	my %counts;
	while (my $row = $Stream->next_row) {
		my $l = $row->length;
		my $c = $row->value(3); # I don't have a score method
		$counts{$c} += $l;
	}
	
	
	# sum bases recorded
	my $sum = sum(values %counts);
	my $hi  = max(keys %counts);
	my $low = min(keys %counts);
	my $ave = mean(keys %counts);
	my $med = median(keys %counts);
	my $std = stddev(keys %counts);
	
	# add to table
	my $i = $Data->add_column($Stream->basename);
	$Data->metadata($i, 'recorded_coverage', "$sum\bp");
	$Data->metadata($i, 'max_coverage', sprintf("%.0f", $hi));
	$Data->metadata($i, 'mean_coverage', sprintf("%.0f", $ave));
	$Data->metadata($i, 'stdev_coverage', sprintf("%.0f", $std));
	$Data->metadata($i, 'median_coverage', sprintf("%.0f", $med));
	$Data->iterate( sub {
		my $row  = shift;
		my $d = $row->value($depth_i); # the depth for the current row
		my $c = 0;
		if ($d == 0) {
			$c = $counts{0}; # count the 0 separately from the other coverage
		}
		else {
			foreach (keys %counts) {
				next if $_ == 0; # do NOT count the 0 coverages here!!!!
				if ($_ >= $d) {
					$c += $counts{$_};
				}
			}
		}
		$row->value($i, sprintf("%.03f", $c / $sum) );
	} );
	print " $sum bases collected from $file\n  min coverage of $low\n  max coverage of $hi\n  mean coverage $ave ± $std\n  median coverage $med\n";
}

my $success = $Data->save($outfile);
print "wrote file $success\n";
