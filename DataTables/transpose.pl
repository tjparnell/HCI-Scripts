#!/usr/bin/perl

use strict;
use Bio::ToolBox::Data;


unless (@ARGV) {
	print "A quick script to transpose a table\nUsage: $0 <input_file>\n";
	exit 0;
}

my $infile = shift (@ARGV);

my $inData = Bio::ToolBox::Data->new(file => $infile);

my $outData = Bio::ToolBox::Data->new(
	feature     => $inData->feature,
);

# transpose the data
for my $r (0 .. $inData->last_row) {
	my $values = $inData->row_values($r);
	$outData->add_column($values);
}

my $outfile = join(q(), $inData->path, $inData->basename, '_transpose', $inData->extension);

# do not include metadata here, not useful
my $success = $outData->write_file(filename => $outfile, simple => 1);

if ($success) {
	print " wrote transposed file $success\n";
}
else {
	warn " could not write file '$outfile'!\n";
}



