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
use Bio::ToolBox;

my $VERSION = 1.6;

unless (@ARGV) {
	print <<END;

A script to summarize the variants in each gene for each sample from 
a VEP annotation table. The output table is columns of samples and 
rows of genes. More than one variant may be found in a particular 
gene, particularly with multi-sample VCFs, and this will summarize 
the table for you.

See the script merge_vcf_samples_with_vep_table.pl to merge a VEP 
annotation (tab-delimited) text file with the sample information 
from the original VCF file. The sample Genotype (GT) should be 
recorded in the file. 

You may wish to filter the table for damaging-only variants prior 
to running this script.

Usage: $0 --input <vep.txt> --output <genes.txt>


END
	exit;
}


my $infile;
my $outfile;

GetOptions( 
	'input=s'   => \$infile, # input table
	'output=s'  => \$outfile, # output table
) or die "bad options!\n";


my $Data = Bio::ToolBox->load_file($infile) or 
	die "unable to open '$infile'!\n";


# Identify columns
my $gene_i = $Data->find_column('gene');
my $symb_i = $Data->find_column('SYMBOL');
my @samples;
{
	# we will look at each value in the first data row and look for something 
	# that looks like a typical genotype 0/1 or ./. something to that effect
	# assume that column is a sample
	my @values = $Data->row_values(1);
	for my $i (0 .. $#values) {
		if ($values[$i] =~ m/[\.01][\/|][\.012]/) {
			push @samples, $i;
		}
	}
}

# Count up
my %gene2count;
my %gene2symbol;
my @gene_order;
my $count = 0;
$Data->iterate( sub {
	my $row = shift;
	
	# get gene information
	my $gene = $row->value($gene_i);
	unless (exists $gene2count{$gene}) {
		$gene2count{$gene} = [ map {0} @samples ];
		$gene2symbol{$gene} = $row->value($symb_i);
		push @gene_order, $gene;
	}
	
	# summarize for gene
	for my $i (0 .. $#samples) {
		my $gt = $row->value($samples[$i]);
		if ($gt =~ m/[01][\/|][12]/) {
			# we have a variant matching 0/1, 1/1, 0/2, or something like that
			$gene2count{$gene}->[$i] += 1;
		}
	}
	$count++;
});


# Write out new table
my $Out = Bio::ToolBox->new_data('Gene', 'Symbol', (map {$Data->name($_)} @samples) );
foreach my $gene (@gene_order) {
	$Out->add_row( [
		$gene,
		$gene2symbol{$gene},
		@{ $gene2count{$gene} }
	] );
}
my $success = $Out->write_file($outfile);

printf " Collected %d genes for %d variants\n", scalar(@gene_order), $count;
print " Wrote file $success\n";





