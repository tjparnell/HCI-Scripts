#!/usr/bin/env perl

use strict;
use Getopt::Long;
use Bio::ToolBox::Data;
use Bio::ToolBox::utility;


### Documentation
unless (@ARGV) {
	print <<END;

A script to calculate ∆∆RPM from four data files. 
Four data files must be generated with get_relative_data (or possibly get_binned_data).
Identical parameters must be employed. Ideally summed normalized read counts are 
collected in windows. A delta of deltas will be calculated from the four input files.

Application: Calculating histone turnover as described in Svensson et al 2015.
    
    (RPKM_T7_2h – RPKM_T7_0h) – (RPKM_HA_2h − RPKM_HA_0h)

Or, in generic terms, using old and new tags at 1st and 2nd time points.
	
	(new_2nd - new_1st) - (old_2nd - old_1st)

A new table will be written.

Usage: $0 --n1 <file> --n2 <file> --o1 <file> --o2 <file> --out <file>

Required: 
  --n1 <file>     New tag, 1st timepoint file
  --n2 <file>     New tag, 2nd timepoint file
  --o1 <file>     Old tag, 1st timepoint file
  --o2 <file>     Old tag, 2nd timepoint file
  --out <file>    The output file.

Optional: 
  --index <range> The range of indices of the columns to calculate.
                  These should be automatically determined by default,
                  unless you have wierd column names


END
	exit;
}


### Options

my $n2_file;
my $o2_file;
my $n1_file;
my $o1_file;
my $outfile;
my $index;

GetOptions( 
	'n2=s'          => \$n2_file, # n2
	'o2=s'          => \$o2_file, # o2
	'n1=s'          => \$n1_file, # n1
	'o1=s'          => \$o1_file, # o1
	'out=s'         => \$outfile, 
	'index=s'       => \$index,
) or die "bad options!\n";

die "missing options! see help" unless 
	($n2_file and $o2_file and $n1_file and $o1_file and $outfile);


### Load input files
print " Loading $n2_file...\n";
my $n2Data = Bio::ToolBox::Data->new( file => $n2_file ) or 
	die " unable to load $n2_file!\n";

print " Loading $o2_file...\n";
my $o2Data = Bio::ToolBox::Data->new( file => $o2_file ) or 
	die " unable to load $o2_file!\n";

print " Loading $n1_file...\n";
my $n1Data = Bio::ToolBox::Data->new( file => $n1_file ) or 
	die " unable to load $n1_file!\n";

print " Loading $o1_file...\n";
my $o1Data = Bio::ToolBox::Data->new( file => $o1_file ) or 
	die " unable to load $o1_file!\n";


### Check files
die " unequal column number!" if 
	$n2Data->number_columns != $o2Data->number_columns or
	$n2Data->number_columns != $n1Data->number_columns or
	$n2Data->number_columns != $o1Data->number_columns;
die " unequal row number!" if 
	$n2Data->last_row != $o2Data->last_row or
	$n2Data->last_row != $n1Data->last_row or
	$n2Data->last_row != $o1Data->last_row;
# I'm assuming the column names are the same too!



### determine indices
my @indices = parse_list($index) if $index;
unless (@indices) {
	for my $i (0 .. $n2Data->number_columns - 1) {
		next if $i eq $n2Data->id_column;
		next if $i eq $n2Data->name_column;
		next if $i eq $n2Data->type_column;
		next if $i eq $n2Data->strand_column;
		next if $i eq $n2Data->chromo_column;
		next if $i eq $n2Data->start_column;
		next if $i eq $n2Data->stop_column;
		push @indices, $i;
	}
	printf " Using indices %s\n", join(',', @indices);
}


### Generate output 
my $outData = $n2Data->duplicate;
foreach my $column (@indices) {
	$outData->metadata($column, 'DeltaDeltaCalculation', join(',', 
		'n2:' . $n2Data->basename,
		'o2:' . $o2Data->basename,
		'n1:' . $n1Data->basename,
		'o1:' . $o1Data->basename,
	) );
}


### Calculate the delta delta
for my $row (1 .. $n2Data->last_row) {
	
	# row values to be stored in the output
	my @values = $n2Data->row_values($row);
	
	# calculate delta for each column index
	for my $column (@indices) {
		# I am assuming these are all valid numbers suitable for subtraction
		# null values should evaluate as zero
		# dot null values will throw an error - I should catch these....
		# store the calculated delta value into the values array
		$values[ $column ] = 
			( $n2Data->value($row, $column) - $n1Data->value($row, $column) ) -
			( $o2Data->value($row, $column) - $o1Data->value($row, $column) );
	}

	# add the new values to the output data array
	$outData->add_row(\@values);
}


### Write output
my $s = $outData->save($outfile);
print " wrote file $s\n";



