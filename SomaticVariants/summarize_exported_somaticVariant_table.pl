#!/usr/bin/env perl

use strict;
use Getopt::Long;
use Bio::ToolBox::Data;

unless (@ARGV) {
	print <<END;

Summarize Somatic variant tables into gene hit count tables. 
Assume tables are from USeq VCFReporter text output. 

If you want to optionally filter based on the maximum normal or delta frequency, 
then run add_annot_table_freqs.pl to put the frequency columns in, then run 
this script. It assumes there are then columns 'NORMAL_AltFrequency' and 
'TUMOR_AltFrequency'.

Usage: $0 <outfile> <file1 file2 ...>
    Options:
     --norm <integer>   maximum normal alt frequency, default 100
     --delta <integer>  minimum normal tumor difference alt frequency, default 0
END
	exit 0;
}

# filtering criteria
my $max_norm_freq = 100;
my $delta_freq = 0;

GetOptions( 
	'norm=i'      => \$max_norm_freq, 
	'delta=i'     => \$delta_freq, 
);

# output file is first from arguments
my $outfile = shift;

# gene hash to store all the results
my %genes;

# walk through all input files
foreach my $file (@ARGV) {
	print "loading file $file...\n";
	my $Data = Bio::ToolBox::Data->new(file => $file) or 
		die "unable to open '$file'\n";
	
	my $norm_freq_i = $Data->find_column('NORMAL_AltFrequency');
	my $tumor_freq_i  = $Data->find_column('TUMOR_AltFrequency');
	my $ensgene_i = $Data->find_column('EnsemblName');
	my $refgene_i = $Data->find_column('RefSeqName');
	my $ensvar_i  = $Data->find_column('EnsemblVarType');
	my $refvar_i  = $Data->find_column('RefSeqVarType');
	
	my $sample;
	if ($Data->basename =~ m/^(.+)_somaticCalls/) {
		$sample = $1;
	}
	elsif ($Data->basename =~ m/^(.+)_filtSomaticCalls/) {
		$sample = $1;
	}
	else {
		warn "unable to identify the sample name! Fix me!!!\n";
	}

	
	my $count = 0;
	$Data->iterate( sub {
		my $row = shift;
		
		# apply filters
		if (defined $norm_freq_i and defined $tumor_freq_i) {
			next if $row->value($norm_freq_i) > $max_norm_freq;
				# maximum normal frequency
			next if ( $row->value($tumor_freq_i) - $row->value($norm_freq_i) ) < $delta_freq;
				# minimum delta frequency
		}
		
		my ($ensgene) = $row->value($ensgene_i) =~ /^(ENSG\d+)/;
		my ($refgene) = $row->value($refgene_i) =~ /^(\w+)/;
		my $variant = $row->value($ensvar_i);
		if ($variant eq 'NA') {
			# to deal with those pesky variants without a ensembl variant type, look 
			# for the refSeq variant type
			unless ($row->value($refvar_i) eq 'NA') {
				$variant = $row->value($refvar_i);
			}
		}
		my $start = $row->start;
		
		unless (exists $genes{$ensgene}) {
			$genes{$ensgene} = {
				'refGene'    => $refgene,
				'count'      => 0,
				'nonsynonymous_SNV' => [],
				'stoploss' => [],
				'stopgain' => [],
				'frameshift' => [],
				'nonframeshift_indel' => [],
				'other' => [],
				'samples' => {},
			};
		}
		
		if ($variant =~ /^nonsynonymous/i) {
			push @{ $genes{$ensgene}{nonsynonymous_SNV} }, "$sample:$start";
		} elsif ($variant =~ /stoploss/) {
			push @{ $genes{$ensgene}{stoploss} }, "$sample:$start";
		} elsif ($variant =~ /stopgain/) {
			push @{ $genes{$ensgene}{stopgain} }, "$sample:$start";
		} elsif ($variant =~ /^nonframeshift/) {
			push @{ $genes{$ensgene}{nonframeshift_indel} }, "$sample:$start";
		} elsif ($variant =~ /^frameshift/) {
			push @{ $genes{$ensgene}{frameshift} }, "$sample:$start";
		} else {
# 			print "other variant '$variant'\n";
			push @{ $genes{$ensgene}{other} }, "$sample:$start";
		}
		$genes{$ensgene}{count}++;
		$genes{$ensgene}{samples}{$sample} += 1;
		
		$count++;
	} );
	
	printf "  Collected $count variants from %s total\n", $Data->last_row;
}

my $Out = Bio::ToolBox::Data->new('columns' => [
	qw(EnsGeneName RefSeqName Total_Number Sample_Number Samples nonsynonymous_SNV_Number Stop_gain_Number 
		Stop_loss_Number Frameshift_Number Nonframeshift_Number Other_Number 
		Nonsynonymous_SNV Stop_gain Stop_loss Frameshift  Nonframeshift Other)
]);

foreach my $gene (keys %genes) {
	
	# get decreasing list of sample IDs and counts
	my @sample_list;
	foreach (
		sort { $genes{$gene}{samples}{$b} <=> $genes{$gene}{samples}{$a} } 
		keys %{ $genes{$gene}{samples} }
	) {
		push @sample_list, "$_:" . $genes{$gene}{samples}{$_};
	}
	
	
	# build an array of values that will become the output row
	my @values = (
		$gene, 
		$genes{$gene}{'refGene'}, 
		$genes{$gene}{count}, 
		scalar(@sample_list),
		join(',', @sample_list),
		scalar @{ $genes{$gene}{nonsynonymous_SNV} }, 
		scalar @{ $genes{$gene}{stopgain} }, scalar @{ $genes{$gene}{stoploss} },
		scalar @{ $genes{$gene}{frameshift} }, 
		scalar @{ $genes{$gene}{nonframeshift_indel} }, 
		scalar @{ $genes{$gene}{other} }, 
	);
	
	# now add sampleIDs for each group
	foreach my $id (
		qw(nonsynonymous_SNV stopgain stoploss frameshift nonframeshift_indel other)
	) {
		push @values, join(',', @{ $genes{$gene}{$id} });
	}
	$Out->add_row(\@values);
}

$Out->sort_data(2, 'd');
my $success = $Out->write_file($outfile);
print "wrote file $success\n";




