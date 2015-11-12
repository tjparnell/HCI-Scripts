#!/usr/bin/env perl

use strict;
use File::Spec;
use Bio::ToolBox::Data;
use constant LOG10 => log(10);

my $bamfile = shift @ARGV;
my (undef, $path, $filename) = File::Spec->splitpath($bamfile);
$filename =~ s/\.bam$//;

my $Data = Bio::ToolBox::Data->new(columns => ["Bin", "$filename\_Count"]);
$Data->database($bamfile);
my $db = $Data->open_database;
my $bam = $db->bam; # open low level bam representation

# get the quality scores from the bam file
my %bin2count;
while (my $a = $bam->read1) {
	$bin2count{ $a->qual }++;
}

# Put counts into table
foreach my $n (0 .. 255) {
	my $c = $bin2count{$n} || 1;
	$Data->add_row( [$n, sprintf("%.2f", log($c)/LOG10) ] );
}

# write file
my $outfile = $bamfile;
$outfile =~ s/\.bam$/_mapq_counts.txt/;
my $success = $Data->save($outfile);
print "wrote file $success\n";
