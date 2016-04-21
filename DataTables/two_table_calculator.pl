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
use Bio::ToolBox::utility;
my $VERSION = 1.1;

### Documentation
unless (@ARGV) {
	print <<END;

A script to compare the cells from two data tables and write a third data table.
The cells can be added, subtracted, multiplied, or divided. 
Alternatively, the minimum or maximum (either as is or absolute values) can be 
taken.

A new table will be written.

Usage: $0 --t1 <file> --t2 <file> --out <file>

Required: 
  --t1 <file>       1st table file
  --t2 <file>       2nd table file
  --method <text>   One of [add|subtract|multiply|divide|min|max|absmin|absmax]
  --out <file>      The output file.

Optional: 
  --index <range> The range of indices of the columns to calculate.
                  These should be automatically determined by default,
                  unless you have wierd column names


END
	exit;
}


### Options

my $t1_file;
my $t2_file;
my $method;
my $outfile;
my $index;

GetOptions( 
	't1=s'          => \$t1_file, # 
	't2=s'          => \$t2_file, # 
	'method=s'      => \$method,
	'out=s'         => \$outfile, 
	'index=s'       => \$index,
) or die "bad options!\n";

die "missing options! see help" unless 
	($t1_file and $t2_file and $method and $outfile);



### check method
my $function;
if ($method eq 'add') {
	$function = sub {
		my $n1 = $_[0] eq '.' ? 0 : $_[0];
		my $n2 = $_[1] eq '.' ? 0 : $_[1];
		return $n1 + $n2;
	}
}
elsif ($method eq 'subtract') {
	$function = sub {
		my $n1 = $_[0] eq '.' ? 0 : $_[0];
		my $n2 = $_[1] eq '.' ? 0 : $_[1];
		return $n1 - $n2;
	}
}
elsif ($method eq 'multiply') {
	$function = sub {
		my $n1 = $_[0] eq '.' ? 0 : $_[0];
		my $n2 = $_[1] eq '.' ? 0 : $_[1];
		return $n1 * $n2;
	}
}
elsif ($method eq 'divide') {
	$function = sub {
		my $n1 = $_[0] eq '.' ? 0 : $_[0];
		my $n2 = $_[1] eq '.' ? 0 : $_[1];
		return '.' if $n2 == 0;
		return $n1 / $n2;
	}
}
elsif ($method eq 'min') {
	$function = sub {
		my $n1 = $_[0] eq '.' ? 0 : $_[0];
		my $n2 = $_[1] eq '.' ? 0 : $_[1];
		return $n1 >= $n2 ? $n2 : $n1;
	}
}
elsif ($method eq 'max') {
	$function = sub {
		my $n1 = $_[0] eq '.' ? 0 : $_[0];
		my $n2 = $_[1] eq '.' ? 0 : $_[1];
		return $n1 >= $n2 ? $n1 : $n2;
	}
}
elsif ($method eq 'absmin') {
	$function = sub {
		my $n1 = $_[0] eq '.' ? 0 : $_[0];
		my $n2 = $_[1] eq '.' ? 0 : $_[1];
		return abs($n1) >= abs($n2) ? $n2 : $n1;
	}
}
elsif ($method eq 'absmax') {
	$function = sub {
		my $n1 = $_[0] eq '.' ? 0 : $_[0];
		my $n2 = $_[1] eq '.' ? 0 : $_[1];
		return abs($n1) >= abs($n2) ? $n1 : $n2;
	}
}
else {
	die "unknown method '$method'! Must use one of [add|subtract|multiply|divide|min|max]\n";
}




### Load input files
print " Loading $t1_file...\n";
my $t1Data = Bio::ToolBox::Data->new( file => $t1_file ) or 
	die " unable to load $t1_file!\n";

print " Loading $t2_file...\n";
my $t2Data = Bio::ToolBox::Data->new( file => $t2_file ) or 
	die " unable to load $t2_file!\n";


### Check files
die " unequal column number!" if $t1Data->number_columns != $t2Data->number_columns;
die " unequal row number!" if $t1Data->last_row != $t2Data->last_row;
# I'm assuming the column names are the same too!



### determine indices
my @indices = parse_list($index) if $index;
unless (@indices) {
	for my $i (0 .. $t1Data->number_columns - 1) {
		next if $i eq $t1Data->id_column;
		next if $i eq $t1Data->name_column;
		next if $i eq $t1Data->type_column;
		next if $i eq $t1Data->strand_column;
		next if $i eq $t1Data->chromo_column;
		next if $i eq $t1Data->start_column;
		next if $i eq $t1Data->stop_column;
		push @indices, $i;
	}
	printf " Using indices %s\n", join(',', @indices);
}


### Generate output 
my $outData = $t1Data->duplicate;
foreach my $column (@indices) {
	$outData->metadata($column, $method . '_Calculation', join(',', 
		't1:' . $t1Data->basename,
		't2:' . $t2Data->basename,
	) );
}


### Calculate the delta delta
for my $row (1 .. $t1Data->last_row) {
	
	# row values to be stored in the output
	my @values = $t1Data->row_values($row);
	
	# calculate delta for each column index
	for my $column (@indices) {
		# I am assuming these are all valid numbers suitable for subtraction
		# null values should evaluate as zero
		# dot null values will throw an error - I should catch these....
		# store the calculated delta value into the values array
		$values[ $column ] = 
			&$function( $t1Data->value($row, $column), $t2Data->value($row, $column) );
	}

	# add the new values to the output data array
	$outData->add_row(\@values);
}


### Write output
my $s = $outData->save($outfile);
print " wrote file $s\n";



