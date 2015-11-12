#!/usr/bin/perl

use strict;
use IO::File;
use File::Spec;

# help
unless (@ARGV) {
	print <<USAGE
A script to parse Picard CollectRnaSeqMetrics output from multiple
files into a single table file.

Two files are written, one for the collected values, the other for the 
coverage, suitable for graphing.

Usage: $0 <output_base> <file1 file2 ...>
USAGE
	;
	exit;
}

# output
my $outfile = shift @ARGV;
die "need input files!\n" unless @ARGV;
$outfile =~ s/\.txt$//i; # remove the extension if present


my @stats_table;
my @coverage_table;

# input
my $count = 0;
foreach my $f (@ARGV) {
	$count++;
	process($count, $f);
}

# average each replicate
average_replicates();

# finish
write_output();
print " merged $count files and wrote $outfile\n";
exit;


sub process {
	my ($count, $file) = @_;
	
	# open input
	my $in = IO::File->new($file, 'r');
	unless ($in) {
		warn "unable to read $file: $!\n";
		return;
	}
	
	# add header
		# use the file name as the column header name
		# first clean it up a bit
	my (undef, undef, $name) = File::Spec->splitpath($file);
	$name =~ s/\.txt$//i;
	$name =~ s/[\._]?metrics?//i;
	if ($count == 1) {
		# first one
		# coverage table
		$coverage_table[0][0] = 'Position';
		$coverage_table[0][1] = $name;
	}
	else {
		$coverage_table[0][$count] = $name;
	}
	
	# walk through file
	my $skip = 1;
	while (my $line = $in->getline) {
		# we want to skip past everything until we reach the coverage table
		
		# skip comment lines
		next if ($line =~ /^#/);
		next if ($line !~ /\w+/); # empty lines
		
		# actually, we want the statistics too, so need to watch for that
		if ($line =~ /^PF_BASES/) {
			# we have the header
			my $line_count = 0;
			if ($count == 1) {
				# first file, so we need this line
				push @stats_table, "NAME\t$line";
			}
			
			# grab the statistics lines
			# internal getline loop
			# bad coding style: Don't do this
			while ($line = $in->getline) {
				last unless $line =~ /\w/; # break out of this loop
				push @stats_table, "$name\t$line";
			}
		}
		
		elsif ($line =~ /^normalized_position/) {
			# we have reached the table header
			
			# internal getline loop
			while ($line = $in->getline) {
				last unless $line =~ /\w/; # break out of this loop
				chomp $line;
				my ($pos, $score) = split /\t/, $line;
				my $row = $pos + 1; 
					# it's odd, they count from 0 to 100, 101 positions!!!?????
					# anyway, use the position as the table index
					# can't use 0 since that will be the header line
				if ($count == 1) {
					$coverage_table[$row][0] = $pos;
					$coverage_table[$row][1] = $score;
				}
				else {
					$coverage_table[$row][$count] = $score;
				}
			}
		}
	}
	$in->close;
}


sub average_replicates {
	return unless $count > 1;
	my $index = $count + 1;
	$coverage_table[0][$index] = 'Mean'; # new header
	foreach my $row (1 .. $#coverage_table) {
		my $sum = 0;
		for my $i (1..$count) {$sum += $coverage_table[$row][$i]}
		$coverage_table[$row][$index] = $sum/$count; # mean
	}
}


sub write_output {
	
	# coverage table
	my $out = IO::File->new("$outfile\_coverage.txt", "w");
	unless ($out) {
		die "cannot write coverage file for $outfile: $!\n";
	}
	foreach my $row (@coverage_table) {
		$out->print(join("\t", @$row) . "\n");
	}
	$out->close;
	
	# statistics table
	my $out = IO::File->new("$outfile\_stats.txt", "w");
	unless ($out) {
		die "cannot write stats file for $outfile: $!\n";
	}
	foreach my $row (@stats_table) {
		$out->print($row);
	}
	$out->close;
}



