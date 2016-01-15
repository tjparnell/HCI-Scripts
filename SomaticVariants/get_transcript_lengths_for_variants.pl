#!/usr/bin/perl

use strict;
use Getopt::Long;
use Statistics::Lite qw(mean max);
use Bio::ToolBox::Data;

unless (@ARGV) {
	print <<END;

A script to merge in transcript and CDS lengths for Variants to assist in 
calculating mutations rates.

Looks for either EnsemblVarDesc or RefSeqVarDesc columns derived from AnnoVar annotation. 
Parses out the Gene and Transcript identifiers from it and looks up the lengths
from the provided tables. Two tables can be provided: one for RefSeq and the other for 
Ensembl annotation. See companion script get_ucsc_transcript_lengths.pl for generating 
the tables.

Inserts two columns with the CDS and Transcript lengths for each reference table 
provided, Ensembl or RefSeq.

Note that multiple transcripts can be listed, and lengths for each will be reported 
in a comma-delimited list. The mean or max can optionally be reported as an additional 
column as well.

Usage:
	get_transcript_lengths_for_variants.pl -i <input.txt> -r <reference_table.txt> -o <output.txt>
	
	-i  a txt input file with the EnsemblVarDesc column
	    try exporting your VCF file with USeq VCFReporter app into a text table
	-r  RefSeq table, from get_ucsc_transcript_lengths.pl
	-e  Ensembl table, from get_ucsc_transcript_lengths.pl
	-x  Include a column for the maximum size
	-m  Include a column for the mean size
	-o  optional output file, default is to overwrite input
	
END
	exit;
}

# options
my $infile;
my $reffile;
my $ensfile;
my $do_max;
my $do_mean;
my $outfile;
GetOptions( 
	'i=s'       => \$infile, # input file
	'r=s'       => \$reffile, # RefSeq table of lengths
	'e=s'       => \$ensfile, # Ensembl table of lengths
	'x!'        => \$do_max, # report max transcript size
	'm!'        => \$do_mean, # report mean transcript size
	'o=s'       => \$outfile, # output vcf
) or die "bad options!\n";



# load input
print "loading $infile...\n";
my $Data = Bio::ToolBox::Data->new(file => $infile) or 
	die "unable to load file $infile!\n";


# reference
print "loading lengths from reference tables....\n";
my $refLengths = load_reference_lengths($reffile);
my $ensLengths = load_reference_lengths($ensfile);
die "unable to load any reference lengths!\n" unless ($refLengths or $ensLengths);


# identify and generate new columns
my ($ensDesc_i, $ensTx_i, $ensCDS_i, $ensName_i, $ensTxMax_i, $ensTxMean_i, 
	$ensCDSMax_i, $ensCDSMean_i) = get_indices('Ensembl');
my ($refDesc_i, $refTx_i, $refCDS_i, $refName_i, $refTxMax_i, $refTxMean_i, 
	$refCDSMax_i, $refCDSMean_i) = get_indices('RefSeq');
 

# look up lengths
print "matching up lengths....\n";
$Data->iterate(\&process_row);


# save
$outfile = $infile unless $outfile;
my $s = $Data->save($outfile);
print "wrote file $outfile\n";



sub load_reference_lengths {
	my $file = shift;
	return unless $file;
	my $Ref = Bio::ToolBox::Data->new(file => $file, stream => 1) or 
		die "unable to load file $file!\n";
	my %id2len;
	# I'm being particularly lazy and hardcoding indices based on 
	# get_ucsc_transcript_lengths
	# transcriptID => hash of tx and cds lengths
	while (my $row = $Ref->next_row) {
		$id2len{ $row->value(1) } = {
			tx  => $row->value(2),
			cds => $row->value(3),
		};
		# also record all values for gene identifier
		$id2len{ $row->value(0) }{tx}  ||= [];
		$id2len{ $row->value(0) }{cds} ||= [];
		push @{$id2len{ $row->value(0) }{tx}}, $row->value(2);
		push @{$id2len{ $row->value(0) }{cds}}, $row->value(3);
	}
	return \%id2len;
}


sub get_indices {
	my $lookup = shift;
	my $Desc_i = $Data->find_column($lookup. 'VarDesc');
	return unless defined $Desc_i;
	my $txi = $Data->add_column($lookup . '_Transcript_length');
	my $cdsi = $Data->add_column($lookup . '_CDS_length');
	
	# optional mean or max columns
	my ($TxMax_i, $TxMean_i, $CDSMax_i, $CDSMean_i);
	if ($do_max) {
		$TxMax_i = $Data->add_column($lookup . '_Transcript_Max_Length');
		$CDSMax_i = $Data->add_column($lookup . '_CDS_Max_Length');
	}
	if ($do_mean) {
		$TxMean_i = $Data->add_column($lookup . '_Transcript_Mean_Length');
		$CDSMean_i = $Data->add_column($lookup . '_CDS_Mean_Length');
	}
	
	# be nice and reorganize the data - this is hard, I should make a new easy method!
	my @additions;
	push @additions, $txi;
	push @additions, $TxMean_i if $TxMean_i;
	push @additions, $TxMax_i if $TxMax_i;
	push @additions, $cdsi;
	push @additions, $CDSMean_i if $CDSMean_i;
	push @additions, $CDSMax_i if $CDSMax_i;
	my @neworder;
	for my $i (0 .. $Desc_i) {push @neworder, $i}
	push @neworder, @additions;
	for my $i ($Desc_i + 1 .. $Data->number_columns - scalar(@additions) -1) {
		push @neworder, $i;
	} 
	$Data->reorder_column(@neworder);
	
	# re-lookup all indices, since everything's been moved around <sigh>
	$Desc_i = $Data->find_column($lookup. 'VarDesc');
	$txi = $Data->find_column($lookup . '_Transcript_length');
	$cdsi = $Data->find_column($lookup . '_CDS_length');
	$TxMax_i = $Data->find_column($lookup . '_Transcript_Max_Length');
	$CDSMax_i = $Data->find_column($lookup . '_CDS_Max_Length');
	$TxMean_i = $Data->find_column($lookup . '_Transcript_Mean_Length');
	$CDSMean_i = $Data->find_column($lookup . '_CDS_Mean_Length');
	my $name_i = $Data->find_column($lookup . 'Name');
	
	# return the new indices for things 
	return ($Desc_i, $txi, $cdsi, $name_i, $TxMax_i, $TxMean_i, 
		$CDSMax_i, $CDSMean_i);
}


sub process_row {
	my $row = shift;
	if ($ensLengths and $ensDesc_i) {
		match_lengths($row, $ensDesc_i, $ensTx_i, $ensCDS_i, $ensName_i, $ensTxMax_i, 
			$ensTxMean_i, $ensCDSMax_i, $ensCDSMean_i, $ensLengths);
	}
	if ($refLengths and $refDesc_i) {
		match_lengths($row, $refDesc_i, $refTx_i, $refCDS_i, $refName_i, $refTxMax_i, 
			$refTxMean_i, $refCDSMax_i, $refCDSMean_i, $refLengths);
	}
}


sub match_lengths {
	my ($row, $Desc_i, $txi, $cdsi, $name_i, $TxMax_i, $TxMean_i, $CDSMax_i, 
		$CDSMean_i, $lengths) = @_;
	my @txlengths;
	my @cdslengths;
	# walk through each comma-delimited variant, usually only 1 but sometimes more
	foreach my $v ( split(",", $row->value($Desc_i)) ) {
		my ($geneid, $txid, $exon, $cdna, $protein) = split(":", $v);
		if (exists $lengths->{$txid}) {
			push @txlengths, $lengths->{$txid}{tx};
			push @cdslengths, $lengths->{$txid}{cds};
		}
		elsif ($geneid eq 'NA' and exists $lengths->{ $row->value($name_i) }) {
			# no specified codon change, try the gene name but only if it exists in table
			$geneid = $row->value($name_i);
			push @txlengths, @{ $lengths->{$geneid}{tx} };
			push @cdslengths, @{ $lengths->{$geneid}{cds} };
		}
		else {
			push @txlengths, 0;
			push @cdslengths, 0;
		}
	}
	
	# record the new values
	$row->value($txi, join(";", @txlengths));
	$row->value($cdsi, join(";", @cdslengths));
	if ($do_max) {
		$row->value($TxMax_i, max(@txlengths));
		$row->value($CDSMax_i, max(@cdslengths));
	}
	if ($do_mean) {
		$row->value($TxMean_i, sprintf("%.0f", mean(@txlengths)));
		$row->value($CDSMean_i, sprintf("%.0f", mean(@cdslengths)));
	}
}




