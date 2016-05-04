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
use Bio::ToolBox::utility qw(
	parse_list
	open_to_read_fh
	open_to_write_fh
);
use Statistics::ANOVA;
use constant LOG10 => log(10);
use File::Temp;
my $VERSION = 1.0;

### Documentation
unless (@ARGV) {
	print <<END;

A script to compare the compare the columns of two or more tables using a 
simple ANOVA test. The values in each corresponding column are compared, 
and a transformed -log10(p-value) is reported for each column test. 

The tables must be identical in terms of the same number of rows and 
columns.Output from get_relative_data.pl is the general intended target.

It's possible to give more than one t2 table, and all of them will be 
compared to t1. A new column is written to the output file for each file 
comparison. 

A new table will be written.

Usage: $0 --t1 <file> --t2 <file> --out <file>

Required: 
  --t1 <file>       1st table file
  --t2 <file>       2nd table file
  --out <file>      The output file.

Optional: 
  --index <range>   The range of indices of the columns to calculate.
                    These should be automatically determined by 
                    default,unless you have wierd column names.
  --transform       Report as -log10(p)
  --qval            convert to q-value, requires R and qvalue lib


END
	exit;
}


### Options

my $t1_file;
my @t2_files;
my $outfile;
my $index;
my $do_transform;
my $do_qval;

GetOptions( 
	't1=s'          => \$t1_file, # 
	't2=s'          => \@t2_files, # more than one can be specified
	'out=s'         => \$outfile, 
	'index=s'       => \$index,
	'transform!'    => \$do_transform,
	'qval!'         => \$do_qval,
) or die "bad options!\n";

die "missing options! see help" unless 
	($t1_file and @t2_files and $outfile);



### Load input files
print " Loading $t1_file...\n";
my $t1Data = Bio::ToolBox::Data->new( file => $t1_file ) or 
	die " unable to load $t1_file!\n";





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
}


### Generate output 
my $outData = Bio::ToolBox::Data->new(
		columns => ['Column','Midpoint',],
);
foreach my $column (@indices) {
	my $midpoint = int(
		# this assumes the column metadata has start and stop
		(($t1Data->metadata($column, 'start') + $t1Data->metadata($column, 'stop') ) / 2) + 0.5
	) or undef; 
	$outData->add_row( [$t1Data->name($column), $midpoint] );
}



### Calculate the p value
my $outcolumn = 2; # set this to the index of the output Data column index
foreach my $t2file (@t2_files) {
	print " Loading $t2file...\n";
	my $t2Data = Bio::ToolBox::Data->new( file => $t2file ) or 
		die " unable to load $t2file!\n";


	# Check files
	die " unequal column number!" if $t1Data->number_columns != $t2Data->number_columns;
	die " unequal row number!" if $t1Data->last_row != $t2Data->last_row;
	# I'm assuming the column names are the same too!
	
	# prepare output column
	$outData->add_column( sprintf("%s_%s_Pval", $t1Data->basename, $t2Data->basename) );
	
	# compare the files
	my $row = 1;
	foreach my $column (@indices) {
		my @t1values = $t1Data->column_values($column);
		my @t2values = $t2Data->column_values($column);
	
		# remove column header
		shift @t1values;
		shift @t2values;
	
		# deal with nulls
		foreach (@t1values) {$_ = 0 if $_ eq '.'}
		foreach (@t2values) {$_ = 0 if $_ eq '.'}
	
		# calculate anova
		my $aov = Statistics::ANOVA->new();
		$aov->load_data( {
			$t1Data->basename => \@t1values,
			$t2Data->basename => \@t2values,
		} );
		$aov->anova(
			independent     => 0,
			parametric      => 1,
			ordinal         => 0,
			tails           => 2,
		);
		my $p = $aov->{_stat}{p_value};
		
		
		# add to output
		$outData->value($row, $outcolumn, $p ? $p : 0 );
		$row++;
	}
	
	$outcolumn++;
}

### Convert to q-value
my $r_exec;
if ($do_qval) {
	# find R
	$r_exec = `which r`;
	chomp $r_exec;
	unless ($r_exec) {
		warn " R executable is not in the path!\n";
		$do_qval = 0;
	}
}
if ($do_qval) {
	print "converting to q-values\n";
	# open a temp file for the p values
	my $p_fh = new File::Temp(
		'UNLINK'   => 0,
		'TEMPLATE' => 'p_valuesXXXXX',
	);
	my $p_file = $p_fh->filename;
	print " writing temporary p value file '$p_file'\n";
	
	# write the p values
	for (my $col = 2; $col < $outData->number_columns; $col++) {
		$outData->iterate( sub {
			my $row = shift;
			$p_fh->printf("%s\n", $row->value($col));
		} );
	}
	$p_fh->close;
	
	# write the R batch command
	my $r_fh = new File::Temp(
		'UNLINK'   => 0,
		'TEMPLATE' => 'r_cmdXXXXX',
	);
	my $r_file = $r_fh->filename;
	my $q_file = $p_file . '.out';
	$r_fh->print("library(qvalue)\n");
	$r_fh->print("p <- scan(\"$p_file\")\n");
	$r_fh->print("qobj <- qvalue(p)\n");
	$r_fh->print("write.qvalue(qobj, file=\"$q_file\")\n");
	$r_fh->close;
	
	# execute
	print " running qvalue package in R...\n";
	system $r_exec, 'CMD', 'BATCH', '--no-save', '--no-restore', $r_file;
	
	# check to see if results are written
	if (-e $q_file and -s $q_file) {
		# success! we have a non-zero output file
		print " found output q file\n";
		
		# add new data columns 
		my @qcolumns;
		for (my $col = 2; $col < $outcolumn; $col++) {
			my $qname = $outData->name($col);
			$qname =~ s/Pval$/Qval/;
			push @qcolumns, $outData->add_column($qname);
			print "  adding q-value column $qname\n";
		}
		
		# open the results file
			# this is a simple text file
			# the file is a simple two column, space delimited
		my $q_fh = open_to_read_fh($q_file);
		my $header = $q_fh->getline;
		# advance the file pointer
		
		# load the q values
		foreach my $col (@qcolumns) {
			print "  loading values for column $col\n";
			$outData->iterate( sub {
				my $row = shift;
				my $line = $q_fh->getline;
				my ($p, $q, $fdr, $pi0) = split /\s/, $line, 4;
				$row->value($col, $q);
			} );
		}
		
		# clean up result files
		$q_fh->close;
		unlink $q_file;
		if (-e "$r_file\.Rout") {
			# R may have written an ouput file
			# we no longer need
			unlink "$r_file\.Rout";
		}
	}
	
	else {
		warn " No qvalue results file from R was written!\n";
		if (-e "$r_file\.Rout") {
			# R may have written an ouput file
			# report this
			warn " there is an R output file printed here\n";
			system 'cat', "$r_file\.Rout";
			# we no longer need
			unlink "$r_file\.Rout";
		}
	}
	
	# clean up
	unlink $p_file;
	unlink $q_file if -e $q_file; 
	unlink $r_file;
	
}


# transform all p and q values
if ($do_transform) {
	$outData->iterate( sub {
		my $row = shift;
		for (my $col = 2; $col < $outData->number_columns; $col++) {
			my $v = $row->value($col);
			unless ($v == 0) {
				$row->value($col, -1 * (log($v) / LOG10) );
			}
		}
	} );
}


### Write output
my $s = $outData->save($outfile);
print " wrote file $s\n";



