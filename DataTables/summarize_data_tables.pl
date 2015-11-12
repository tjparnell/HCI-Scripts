#!/usr/bin/perl

# documentation at end of file

use strict;
use Getopt::Long;
use Pod::Usage;
use Statistics::Lite qw(sum mean median min max stddevp);
use Bio::ToolBox::Data;
use Bio::ToolBox::utility qw(ask_user_for_index format_with_commas parse_list);
my $VERSION =  '1';

print "\n This script will summarize two or more data files\n\n";


### Quick help
unless (@ARGV) { 
	# when no command line options are present
	# print SYNOPSIS
	pod2usage( {
		'-verbose' => 0, 
		'-exitval' => 1,
	} );
}


### Get command line options and initialize values
my (
	$outfile,
	$index,
	$method,
	$gz,
	$help,
	$print_version,
);

# Command line options
GetOptions( 
	'out=s'     => \$outfile, # specify the input data file
	'index=s'   => \$index, # specify the index
	'method=s'  => \$method, # specify the method
	'gz!'       => \$gz, # compress output files
	'help'      => \$help, # request help
	'version'   => \$print_version, # print the version
) or die " unrecognized option(s)!! please refer to the help documentation\n\n";

# Print help
if ($help) {
	# print entire POD
	pod2usage( {
		'-verbose' => 2,
		'-exitval' => 1,
	} );
}

# Print version
if ($print_version) {
	print " Biotoolbox script summarize_data_tables.pl, version $VERSION\n\n";
	exit;
}


### Defaults
unless ($outfile) {
	die "must provide an output filename!\n";
}
my @indices;
if ($index) {
	@indices = parse_list($index);
}
if ($method) {
	unless ($method =~ /mean|median|sum|min|max|stdev/) {
		die "unrecognized method '$method'\n";
	}
}
else {
	$method = 'mean';
}
my $methodsub = 
	$method eq 'mean'   ? \&mean   :
	$method eq 'median' ? \&median :
	$method eq 'min'    ? \&min    :
	$method eq 'max'    ? \&max    :
	$method eq 'sum'    ? \&sum    :
	$method eq 'stdev'  ? \&stddevp : undef;


### Output
my $Out = Bio::ToolBox::Data->new();
$Out->add_column('File'); # index 0
$Out->add_column('Method'); # index 1



### Input
my @names; # a global array of the column names so that we can warn user if files differ
while (@ARGV) {
	
	# load input file
	my $file = shift @ARGV;
	print "  summarizing file $file...\n";
	my $Data = Bio::ToolBox::Data->new(file => $file);
	unless ($Data) {
		warn "unable to load file $file!\n";
		next;
	}
	
	# get indices
	unless (@indices) {
		@indices = ask_user_for_index($Data, "Enter one or more indices to summarize  ");
		die "must provide valid indices!\n" unless (@indices);
	}
	
	# add columns to output
	if ($Out->number_columns == 2) {
		foreach (@indices) {
			$Out->add_column( $Data->name($_) );
			$names[$_] = $Data->name($_); # make the names array match the input table
		}
	}
	
	# start building row data
	my @row_data = ($Data->basename, $method);
	
	# collect the data 
	foreach my $i (@indices) {
		unless ($Data->name($i)) {
			# uh oh, these files don't match!? print a friendly warning and record null
			print "   Warning: column $i doesn't exist!\n";
			push @row_data, '.';
			next;
		}
		if ($Data->name($i) ne $names[$i]) {
			# print a friendly warning, but go ahead and proceed
			printf "   Warning: column $i name %s doesn't match %s\n", $Data->name($i), 
				$names[$i];
		}
		my @numbers;
		foreach ($Data->column_values($i)) {
			push @numbers, $_ if /^\-?\d+\.?\d*$/; # looks like an integer or float, no exponent
		}
		push @row_data, &$methodsub(@numbers);
	}
	
	# add to output data
	$Out->add_row(\@row_data);
}

# write output file
my $success = $Out->write_file(filename => $outfile, gz => $gz);
print "wrote file $success\n";


__END__

=head1 NAME

summarize_data_tables.pl

A script to join two or more data files and concatenate rows.

=head1 SYNOPSIS

summarize_data_tables.pl [--options] --out <filename> <file1> <file2> ...
  
  Options:
  --out <filename>
  --index <integer,range>
  --method [sum mean median min max stddevp]
  --gz
  --version
  --help

=head1 OPTIONS

The command line flags and descriptions:

=over 4

=item --out <filename>

Provide the name of the output file. 

=item --gz

Indicate whether the output files should be compressed 
with gzip. 

=item --version

Print the version number.

=item --help

Display the POD documentation

=back

=head1 DESCRIPTION

This program will summarize two or or more data files....

=head1 AUTHOR

 Timothy J. Parnell, PhD
 Huntsman Cancer Institute
 University of Utah
 Salt Lake City, UT, 84112

This package is free software; you can redistribute it and/or modify
it under the terms of the Artistic License 2.0.  

