#!/usr/bin/env perl

use strict;
use Getopt::Long;
use Bio::ToolBox::Data;
use Bio::ToolBox::utility qw(parse_list ask_user_for_index format_with_commas);

my $VERSION = 3;

=head1 CHANGES
version 2

    * Removed filtering based on normal or delta frequency
      This should be done at the VCF level. See filter_vcf_by_FA script.
    * Removed both Ensembl and RefSeq gene id reporting
    * Now pick explicitly either Ensembl or RefSeq gene id reporting
    * No longer assumes one sample per file. Multiple samples per file are 
      taken. Even multiple samples in only one file. Indices of sample columns 
      must be provided or selected interactively from list
    * Allows counting of one or more sample classes based on a text substring 
      in the sample identifiers.

version 3
	* Add support for VEP annotation tables

=cut

unless (@ARGV) {
	print <<END;

Summarize Somatic variant tables into gene hit count tables.  

Tables from USeq VCFReporter text output are supported.

Tables from VEP should also be supported. Filter on the Pick'ed 
consequence if more than one is reported. It is best to use tables 
previously filtered for damaging consequences.

The output table lists all of the genes with 1 or more variant, and 
lists how many and which samples hit each gene. Hits are categorized by 
their type of effect, e.g. nonsynonymous, stop gain, etc. This is useful 
when you have lots of samples and you want to identify the top genes that 
are repeatedly hit.

Tables by default use the Ensembl identifier, or it can optionally use 
the RefSeq gene identifier. The columns EnsemblName and EnsemblVarType, 
or RefSeqName and RefSeqVarType, are assumed to be present for the USeq 
VCFReporter tables. For VEP tables, it uses the SYMBOL and Consequence 
columns.

If there are multiple sample classes and you want individual counts for 
each sample class, provide one or more text substrings to identify and 
count the classes. For example, you have sample identifiers are appended 
with T and AD for tumor and adenoma. Provide T and AD to separate --count 
options to count these. A total count is always provided.

Usage: summarize_exported_somaticVariant_table.pl -o <outfile> <file1 file2 ...>
    Options:
     --out <outfile>    Provide an output name
     --sample <index>   provide the 0-based index, index range, or 
                        comma-delimited list of indices for the columns 
                        with the sample information. If not provided, 
                        they are interactively selected from a list.
     --refseq           use the RefSeq identifier instead of Ensembl
     --count <string>   provide a text substring to identify classes 
END
	exit 0;
}




### OPTIONS
my $outfile = 'summarized_variants.txt';
my $do_refseq = 0;
my $sample_string;
my @count_tags;

GetOptions( 
	'out=s'       => \$outfile,
	'sample=s'    => \$sample_string,
	'refseq!'     => \$do_refseq,
	'count=s'     => \@count_tags,
) or die "something's wrong with your options! check your command\n";

# gene hash to store all the results
# sometimes there are simple genes in the gene id column, sometimes not
# sometimes multiple genes are listed, either because they overlap or variant is intergenic
my %genes;
my %nongenes;



### INPUT FILES
# walk through all input files
die "no input files!\n" unless (@ARGV);
foreach my $file (@ARGV) {
	print "loading file $file...\n";
	my $Data = Bio::ToolBox::Data->new(file => $file) or 
		die "unable to open '$file'\n";
	
	# identify columns
	my $ensgene_i = $Data->find_column('EnsemblName') || 
					$Data->find_column('SYMBOL') || 
					$Data->find_column('Gene') || undef;
	my $refgene_i = $Data->find_column('RefSeqName') ||
					$Data->find_column('SYMBOL') || 
					$Data->find_column('Gene') || undef;
	my $ensvar_i  = $Data->find_column('EnsemblVarType') ||
					$Data->find_column('Consequence');
	my $refvar_i  = $Data->find_column('RefSeqVarType') ||
					$Data->find_column('Consequence');
	my $location_i= $Data->find_column('Location');
	
	# check
	if ($do_refseq and not defined $refgene_i and not defined $refvar_i) {
		die "no RefSeqName and/or RefSeqVarType columns identified in '$file'!\n";
	}
	elsif (not defined $ensgene_i and not defined $ensvar_i) {
		die "no EnsemblName and/or EnsemblVarType columns identified in '$file'!\n";
	}
	
	# identify sample(s)
	my @samples;
	if ($sample_string) {
		@samples = parse_list($sample_string);
		foreach (@samples) {
			die "invalid index $_ for file $file!!!\n" unless ($Data->name($_));
		}
	}
	else {
		@samples = ask_user_for_index($Data, 'Enter the column index for the samples  ');
	}

	# parse the samples
	my $count = 0;
	$Data->iterate( sub {
		my $row = shift;
		
		# get gene identifiers
		my $gene = $do_refseq ? $row->value($refgene_i) : $row->value($ensgene_i);
		my $variant = $do_refseq ? $row->value($refvar_i) : $row->value($ensvar_i);
		
		# pick the appropriate gene hash to store the data in
		my $hash;
		if ($gene =~ /[,\(\)\:]/) {
			# non-gene identifier characters imply something complicated going on
			$hash = \%nongenes;
		}
		else {
			# normal happy gene identifier
			$hash = \%genes;
		}
		
		# build the gene hash for storing
		unless (exists $hash->{$gene}) {
			$hash->{$gene} = {
				'count'      => 0,
				'count_tags' => {},
				'nonsynonymous_SNV' => [],
				'stoploss' => [],
				'stopgain' => [],
				'frameshift' => [],
				'nonframeshift_indel' => [],
				'other' => [],
				'samples' => {},
			};
		}
	
		# add each sample information to the gene hash
		my $start = $row->start;
		if (not $start and $location_i) {
			(undef, $start) = split(':', $row->value($location_i));
		}
		foreach my $i (@samples) {
			my $sampleid = $Data->name($i); # column name is the sample ID
			my ($allele, $gt, undef) = split(/:/, $row->value($i)); 
				# assuming the VCFReporter format here allele:genotype:ref,alt
		
			# check the sample, only count if it is variant
			next if $row->value($i) eq 'NA';
			next if $gt eq '0';
			next if $gt eq '0/0';
		
			# determine variant category
			if ($variant =~ /nonsynonymous/i) {
				push @{ $hash->{$gene}{nonsynonymous_SNV} }, "$sampleid:$start";
			} elsif ($variant =~ /missense.?variant/i) {
				push @{ $hash->{$gene}{nonsynonymous_SNV} }, "$sampleid:$start";
			} elsif ($variant =~ /stop.?loss/) {
				push @{ $hash->{$gene}{stoploss} }, "$sampleid:$start";
			} elsif ($variant =~ /stop.?gain/) {
				push @{ $hash->{$gene}{stopgain} }, "$sampleid:$start";
			} elsif ($variant =~ /non.?frame.?shift/) {
				push @{ $hash->{$gene}{nonframeshift_indel} }, "$sampleid:$start";
			} elsif ($variant =~ /frame.?shift/) {
				push @{ $hash->{$gene}{frameshift} }, "$sampleid:$start";
			} else {
				push @{ $hash->{$gene}{other} }, "$sampleid:$start";
			}
		
			# update counts
			$hash->{$gene}{count}++;
			$hash->{$gene}{samples}{$sampleid} += 1;
			foreach my $tag (@count_tags) {
				if ($sampleid =~ /$tag/) {
					$hash->{$gene}{count_tags}{$tag} += 1;
				}
			}
		
			$count++; # total count for this file
		}
	} );
	
	printf "  Collected %s variants from %s total\n", format_with_commas($count), 
		format_with_commas($Data->last_row);
}



### OUTPUT DATA table
# prepare output file
my @headers = qw(GeneName Total_Number); 
push @headers, map {$_ . '_count'} @count_tags if @count_tags;
push @headers, qw(Sample_Number Samples nonsynonymous_SNV_Number Stop_gain_Number 
		Stop_loss_Number Frameshift_Number Nonframeshift_Number Other_Number 
		Nonsynonymous_SNV Stop_gain Stop_loss Frameshift  Nonframeshift Other);
my $Out = Bio::ToolBox::Data->new('columns' => \@headers);

# start putting the collected gene data into the output data array
# sort first by decreasing total count, then asciibetically by name
# do normal genes first, then the oddball guys
printf "Sorting %s genes\n", format_with_commas(scalar keys %genes);
foreach my $gene (
	map {$_->[1]}
	sort {$b->[0] <=> $a->[0] or $a->[1] cmp $b->[1]}
	map { [$genes{$_}{count}, $_] }
	keys %genes
) {
	add_gene_to_output(\%genes, $gene);
}
printf "Sorting %s other multiple or non-genes\n", format_with_commas(scalar keys %nongenes);
foreach my $gene (
	map {$_->[1]}
	sort {$b->[0] <=> $a->[0] or $a->[1] cmp $b->[1]}
	map { [$nongenes{$_}{count}, $_] }
	keys %nongenes
) {
	add_gene_to_output(\%nongenes, $gene);
}



### WRITE OUTPUT
printf "Summarized %s %s genes\n", format_with_commas($Out->last_row), 
	$do_refseq ? 'RefSeq' : 'Ensembl';

my $success = $Out->write_file($outfile);
print "wrote file $success\n";




### SUBROUTINES

sub add_gene_to_output {
	my ($hash, $gene) = @_;
	
	# get decreasing list of sample IDs and counts
	my @sample_list;
	foreach (
		sort { $hash->{$gene}{samples}{$b} <=> $hash->{$gene}{samples}{$a} } 
		keys %{ $hash->{$gene}{samples} }
	) {
		push @sample_list, "$_:" . $hash->{$gene}{samples}{$_};
	}
	my @tag_counts;
	foreach (@count_tags) {
		push @tag_counts, $hash->{$gene}{count_tags}{$_} || 0;
	}
	
	# build an array of values that will become the output row
	my @values = (
		$gene, 
		$hash->{$gene}{count},
	);
	push @values, @tag_counts if (@count_tags);
	push @values, scalar(@sample_list);
	push @values, join(',', @sample_list);
	push @values, scalar @{ $hash->{$gene}{nonsynonymous_SNV} };
	push @values, scalar @{ $hash->{$gene}{stopgain} };
	push @values, scalar @{ $hash->{$gene}{stoploss} };
	push @values, scalar @{ $hash->{$gene}{frameshift} };
	push @values, scalar @{ $hash->{$gene}{nonframeshift_indel} };
	push @values, scalar @{ $hash->{$gene}{other} };
	
	# now add sampleIDs for each group
	foreach my $id (
		qw(nonsynonymous_SNV stopgain stoploss frameshift nonframeshift_indel other)
	) {
		push @values, join(',', @{ $hash->{$gene}{$id} });
	}
	
	# finally add these values to the output array
	$Out->add_row(\@values);
}

