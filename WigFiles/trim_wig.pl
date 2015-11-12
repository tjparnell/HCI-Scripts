#!/bin/env perl

use warnings;
use strict;

# script to cut low scores from a wig file

$^I = '.bak'; # set the backup file extension
my $total = 0;
my $kept = 0;

while (<>) {
	$total++;
	if (/^variableStep/) {
		print;
		$kept++;
		next;
	}
	my ($p, $s) = split /\s/;
	chomp $s;
	next if ($s < 0.01);
	print;
	$kept++;
}

print " kept $kept lines out of $total lines\n";
