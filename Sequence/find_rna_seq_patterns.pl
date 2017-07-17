#!/usr/bin/perl

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
use Bio::ToolBox::db_helper 1.51 qw(get_genomic_sequence);
use Bio::ToolBox::parser::gff;
use Bio::ToolBox::GeneTools qw(get_exons get_transcripts get_cds get_utrs);

my $VERSION = 1;

unless (@ARGV) {
	print <<END;
 A script to search for sequence patterns or motifs in transcript annotation.
 At this time it's only using a simple, unweighted, regular expression.
 It will report the gene name, transcript id, genomic coordinates, matched 
 sequence, and where in the transcript it was found.
 
 Usage:
 $0 -i <gtf> -o <outfile> -f fasta -p <regex> 
 
 Options
 --in       Input annotation file, GTF or GFF3
 --out      Ouput file
 --fasta    genomic fasta file
 --pattern  A regex pattern to search

END
	exit;
}


### Get Options
my ($infile, $outfile, $pattern, $fasta);
GetOptions( 
	'in=s'       => \$infile, # the input GTF file path
	'out=s'      => \$outfile, # name of output file 
	'pattern=s'  => \$pattern, # the regex to search
	'fasta=s'    => \$fasta, # a searchable fasta file
) or die " unrecognized option(s)!! please refer to the help documentation\n\n";



### Prepare output structure
my $Data = Bio::ToolBox::Data->new(
	columns => [ qw(GeneName GeneID TranscriptID Biotype Subfeature TranscriptPosition Chromosome Start Stop Strand Sequence) ],
);

### Open fasta
my $db = $Data->open_database($fasta) or 
	die "unable to open fasta file! $!\n";


### open and parse GFF file
my $gff = Bio::ToolBox::parser::gff->new($infile) or 
	die "unable to open '$infile' $!\n";
$gff->parse_file;


### Iterate through genes
print " Searching for sequence pattern $pattern....\n";
my %found;
while (my $gene = $gff->next_top_feature) {
	
	# Work on each transcript
	foreach my $transcript (get_transcripts($gene)) {
		
		# assemble sequence and collect coordinates
		my $seq;
		my @gcoords; # genomic coordinates
		my @tcoords; # relative transcript coordinates
		my $length = 0;
		foreach my $e (get_exons($transcript)) {
			$seq .= get_genomic_sequence($db, $e->seq_id, $e->start, $e->end);
			push @gcoords, [$e->start, $e->end];
			push @tcoords, [$length + 1, $length + $e->length];
			$length += $e->length;
		}
		my @ccoords; # cds genomic coordinates
		foreach my $c (get_cds($transcript)) {
			push @ccoords, [$c->start, $c->end];
		}
		my @ucoords; # utr genomic coordinates
		foreach my $u (get_utrs($transcript)) {
			push @ucoords, [$u->start, $u->end, $u->primary_tag];
		}
		
		# search for the pattern
		# we are using funky regex variables @- and @+ that record the positions of match
		my @results; # arrays of [start, stop, sequence]
		if ($transcript->strand < 0) {
			# reverse strand
			$seq =~ tr/gatcGATC/ctagCTAG/;
			$seq = reverse $seq;
			while ($seq =~ m/$pattern/g) {
				push @results, [ $length - $+[0] + 1, $length - $-[0], 
					substr($seq, $-[0], $+[0] - $-[0]) ];
			}
		}
		else {
			# forward strand
			while ($seq =~ m/$pattern/g) {
				push @results, [ $-[0] + 1, $+[0], substr($seq, $-[0], $+[0] - $-[0]) ];
			}
		}
		
		# record the pattern results
		foreach my $r (@results) {
			my ($start, $end, $matchseq) = @$r;
			
			# convert to genomic coordinates
			my ($gstart, $gstop);
			foreach my $i (0 .. $#tcoords) {
				if ($start >= $tcoords[$i][0] and $start <= $tcoords[$i][1]) {
					$gstart = $gcoords[$i][0] + ($start - $tcoords[$i][0]);
				}
				if ($end >= $tcoords[$i][0] and $end <= $tcoords[$i][1]) {
					$gstop = $gcoords[$i][0] + ($end - $tcoords[$i][0]);
				}
				last if ($gstart and $gstop);
			}
			
			# identify the subfeature
			my $type = 'exon'; # if all else fails, it's an exon
			foreach my $c (@ccoords) {
				if ($gstart >= $c->[0] and $gstart <= $c->[1]) {
					$type = 'cds';
					last;
				}
			}
			foreach my $u (@ucoords) {
				if ($gstart >= $u->[0] and $gstart <= $u->[1]) {
					$type = $u->[2];
					last;
				}
			}
			
			# fractional position in transcript
			my $fraction = ($gene->strand > 0) ? $start/$length : ($length-$start)/$length;
			
			# record
			# Gene Transcript Biotype Subfeature Offset_Fraction Chromosome Start Stop Strand Sequence
			$Data->add_row( [
				$gene->display_name,
				$gene->primary_id,
				$transcript->primary_id,
				$transcript->primary_tag,
				$type,
				sprintf("%.3f", $fraction),
				$gene->seq_id,
				$gstart,
				$gstop,
				$gene->strand > 0 ? '+' : '-',
				$matchseq,
			] );
			
			# count of found matches
			$found{$gene->primary_id} += 1;
		}
	}
}
printf " Found %d sequences in %d genes\n", $Data->last_row, scalar(keys %found);



### Sort the table nicely
print " Sorting found sequences\n";
$Data->gsort_data;


### Write data
my $s = $Data->save($outfile);
printf " Wrote file $s\n";


