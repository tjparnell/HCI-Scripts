#!/usr/bin/perl

use strict;
use Getopt::Long;
use Bio::ToolBox::Data;
use Bio::ToolBox::parser::ucsc; # change to gff as necessary

unless (@ARGV) {
	print <<END;

A script to print out the transcript and CDS lengths of transcripts.

Currently assuming a UCSC-formatted gene table, but could also use GFF3 or GTF.

Prints out gene identifier, transcript identifier, transcript length, CDS length.
Note that non-coding genes will have no CDS length!

Usage:
	$0 -i <ucsc> -o <output.txt>

END
	exit;
}

# options
my $infile;
my $outfile;
GetOptions( 
	'i=s'       => \$infile, # input file
	'o=s'       => \$outfile, # output vcf
) or die "bad options!\n";



# parse the table into seqfeatures
print "parsing $infile....\n";
my $ucsc = Bio::ToolBox::parser::ucsc->new(
        file      => $infile,
        do_cds    => 1, # cds is by default off for ucsc files, but we want this
) or die "unable to parse $infile!\n";



# collect the data
my $Data = Bio::ToolBox::Data->new(
	'columns'  => [qw(GeneID TranscriptID Transcript_Length CDS_Length)],
);
while (my $gene = $ucsc->next_top_feature) {
	my $geneName = $gene->display_name;
	foreach my $t (get_transcripts($gene)) {
		my $txlength = get_transcript_length($t);
		my $cdslength = get_transcript_cds_length($t);
		$Data->add_row( [$geneName, $t->display_name, $txlength, $cdslength] );
	}
}


# finish
printf "collected %d transcripts from $infile\n", $Data->last_row;
my $s = $Data->save($outfile);
print "wrote file $s\n" if $s;




#### methods from github branch seqfeature not yet incorporated into BioToolBox

# currently package Bio::ToolBox::Gene::Utility, but subject to change....
# just a subset of the methods I think I need....
# ideally I would be calling just one db parser

sub get_exons {
	
	# initialize
	my $transcript = shift;
	my @exons;
	my @cdss;
	my @transcripts;
	
	# go through the subfeatures
	foreach my $subfeat ($transcript->get_SeqFeatures) {
		my $type = $subfeat->primary_tag;
		if ($type =~ /exon/i) {
			push @exons, $subfeat;
		}
		elsif ($type =~ /cds|utr|untranslated/i) {
			push @cdss, $subfeat;
		}
		elsif ($type =~ /rna|transcript/i) {
			push @transcripts;
		}
	}
	
	# check which array we'll use
	# prefer to use actual exon subfeatures, but those may not be defined
	my @list;
	if (@exons) {
		@list = @exons;
	}
	elsif (@cdss) {
		@list = @cdss;
	}
	elsif (@transcripts) {
		foreach my $t (@transcripts) {
			# there are possibly duplicates in here if there are alternate transcripts
			# should we remove them?
			push @list, collect_exons($t);
		}
	}
	else {
		# nothing found!
		return;
	}
	
	# return sorted list by start position
	return  map { $_->[0] }
			sort { $a->[1] <=> $b->[1] }
			map { [$_, $_->start] } 
			@list;
}

sub get_transcript_length {
	my $transcript = shift;
	my $total = 0;
	foreach my $e (get_exons($transcript)) {
		$total += $e->length;
	}
	return $total;
}

sub get_transcript_cds_length {
	my $transcript = shift;
	my $total = 0;
	foreach my $subf ($transcript->get_SeqFeatures) {
		next unless $subf->primary_tag eq 'CDS';
		$total += $subf->length;
	}
	return $total;
}

sub get_transcripts {
	my $gene = shift;
	my @transcripts;
	foreach my $subf ($gene->get_SeqFeatures) {
		push @transcripts, $subf if $subf->primary_tag =~ /rna|transcript/i;
	}
	return map { $_->[0] }
		sort { $a->[1] <=> $b->[1] }
		map { [$_, $_->start] } 
		@transcripts;
}


