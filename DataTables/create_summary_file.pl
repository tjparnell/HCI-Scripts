#!/usr/bin/perl

use strict;
use Bio::ToolBox::Data;

print "Creates a summary file from a data file using the builtin Data summary_file method\n";
unless (@ARGV) {
	print "Usage: $0 <file1> <file2> ...\n";
	exit;
}

foreach my $file (@ARGV) {
	my $Data = Bio::ToolBox::Data->new(file => $file);
	unless ($Data) {
		warn "can't process $file! $!\n";
		next;
	}
	my $sumfile = $Data->summary_file();
	print "wrote summary file $sumfile\n";
}


