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
use Bio::ToolBox::Data '1.40';
my $VERSION = 1.5;

unless (@ARGV) {
	print <<END;

A script to fix and standardize sample attributes in somatic VCF files. 
Should properly handle somatic VCF files from MuTect, SomaticSniper, 
Strelka, VarScan, and scalpel.

It performs the following functions:
  - adds missing alternate fraction (FA or FREQ) attribute (Strelka and 
    SomaticSniper). Changes FREQ to FA (VarScan)
  - adds missing genotype GT attribute (Strelka)
  - adds ref,alt allele depth (AD) attribute (Strelka, VarScan)
  - optionally discard excess sample attributes, leaving only GT:AD:FA
    otherwise merging VCFs from different callers with different sample 
    attributes into a single VCF often leads to malformed VCF sample records

Usage: $0 -i <input.vcf> 
Options: 
  -i <input.vcf>        input vcf file, may be gz compressed
  -o <output.vcf>       default is to overwrite input!!!!
  -t <NAME>             Tumor sample name, default TUMOR
  -n <NAME>             Normal sample name, default NORMAL
  -m                    flag to minimize to just three sample attributes 
                        (GT:AD:FA) and remove all INFO fields
END
	exit;
}

my $file;
my $outfile;
my $nName = 'NORMAL';
my $tName = 'TUMOR';
my $minimize;

GetOptions( 
	'i=s'       => \$file, # input vcf
	'o=s'       => \$outfile, # output vcf
	'n=s'       => \$nName, # normal name
	't=s'       => \$tName, # tumor name
	'm!'        => \$minimize, # throw away extra attributes
) or die "bad options!\n";

# unnecessary INFO fields to remove when minimizing
# this list may need to be updated with keys from other callers
# not updating will leave behind those INFO fields, not a huge problem


# input file
my $Data = Bio::ToolBox::Data->new(in => $file) or 
	die "unable to open $file!\n";
die "file is not VCF!!!\n" unless $Data->vcf;

# check names
my $nNameIndex = $Data->find_column($nName);
my $tNameIndex = $Data->find_column($tName);
unless ($nNameIndex and $tNameIndex) {
	print "Unable to identify Normal sample '$nName' column!\n" unless $nNameIndex;
	print "Unable to identify Tumor sample '$tName' column!\n" unless $tNameIndex;
	print "Please check your VCF file and provide the correct sample names.\n";
	exit;
}


# output file
$outfile ||= $file;
$outfile =~ s/\.gz$//i;
$outfile = $outfile . '.vcf' unless $outfile =~ /\.vcf$/i;



# process
my $change_AD_comment = 0; # remember to add the AD header if necessary
my $iterator = $Data->row_stream;
while (my $line = $iterator->next_row) {
	
	my $att = $line->vcf_attributes;
	
	# Strelka data tier and variant type, used in both FA and depth calculations
	my ($tier, $indel);
	
	### check for FA or FREQ
	if (exists $att->{$nNameIndex}{FREQ}) {
		# VarScan somatic attribute, convert to MuTect style FA for consistency
		my $nFA = $att->{$nNameIndex}{FREQ};
		$nFA =~ s/\$$//;
		$nFA /= 100;
		$att->{$nNameIndex}{FA} = $nFA;
		my $tFA = $att->{$tNameIndex}{FREQ};
		$tFA =~ s/\$$//;
		$tFA /= 100;
		$att->{$tNameIndex}{FA} = $tFA;
	}
	unless (exists $att->{$nNameIndex}{FA}) {
		# No frequency of alternate
		
		# calculate FREQ attribute
		my ($normal_FA, $tumor_FA);
		if (exists $att->{$nNameIndex}{DP4}) {
			# DP4 = ref-forward bases, ref-reverse, alt-forward and alt-reverse bases
			my ($rf, $rr, $af, $ar) = split ',', $att->{$nNameIndex}{DP4}, 4;
			$normal_FA = sprintf "%0.2f", ($af + $ar) / ($rf + $rr + $af + $ar);
			($rf, $rr, $af, $ar) = split ',', $att->{$tNameIndex}{DP4}, 4;
			$tumor_FA = sprintf "%0.2f", ($af + $ar) / ($rf + $rr + $af + $ar);
		} 
		
		elsif (exists $att->{$nNameIndex}{TIR} or exists $att->{$nNameIndex}{TAR}) {
			# Strelka attributes
			
			# determine data tier and type of indel
			if (exists $att->{INFO}{TQSI}) {
				# strelka indel
				$tier = $att->{INFO}{TQSI};
				$indel = 1;
			}
			elsif (exists $att->{INFO}{TQSS}) {
				# strelka somatic
				$tier = $att->{INFO}{TQSS};
				$indel = 0;
			}
			unless ($tier) {
				# just go with the first data tier 
				warn sprintf "no tier data for data line %s:%d\n", 
					$line->seq_id, $line->start;
				$tier = 1;
			}
			
			# get the supporting number of alternate reads based on data tier and 
			# somatic variant type
			my ($nt1, $nt2) = split ',', # normal
								$indel ? $att->{$nNameIndex}{TIR} : $att->{$nNameIndex}{TAR}, 
								2;
			my ($tt1, $tt2) = split ',', # tumor
								$indel ? $att->{$tNameIndex}{TIR} : $att->{$tNameIndex}{TAR}, 
								2;
			
			# calculate FREQ based on which data tier was used
			# using the total depth DP as the denominator
			if ($tier == 1) {
				$normal_FA = $att->{$nNameIndex}{DP} ? 
					sprintf "%0.2f", $nt1 / $att->{$nNameIndex}{DP} : 0;
				$tumor_FA  = $att->{$tNameIndex}{DP} ? 
					sprintf "%0.2f", $tt1 / $att->{$tNameIndex}{DP} : 0;
			} elsif ($tier == 2) {
				$normal_FA = $att->{$nNameIndex}{DP2} ? 
					sprintf "%0.2f", $nt1 / $att->{$nNameIndex}{DP2} : 0;
				$tumor_FA  = $att->{$tNameIndex}{DP2} ? 
					sprintf "%0.2f", $tt1 / $att->{$tNameIndex}{DP2} : 0;
			} else {
				warn sprintf "unknown tier '$tier' for data line %s:%d\n",
					$line->seq_id, $line->start;
				$normal_FA = 0;
				$tumor_FA  = 0;
			}
		} 
		
		elsif (exists $att->{$nNameIndex}{AD}) {
			# AD = ref bases, alt bases
			my ($ref, $alt) = split ',', $att->{$nNameIndex}{AD}, 2;
			my $total = $ref + $alt;
			$normal_FA = $total == 0 ? '0.00' : sprintf "%0.2f", $alt / $total;
			($ref, $alt) = split ',', $att->{$tNameIndex}{AD}, 2;
			$total = $ref + $alt;
			$tumor_FA = $total == 0 ? '0.00' : sprintf "%0.2f", $alt / $total;
		} 
		
		else {
			warn sprintf "missing read depth identifiers for data line %s:%d\n",
				$line->seq_id, $line->start;
			$normal_FA = 0;
			$tumor_FA  = 0;
		}
		
		# push the new FA attribute
		$att->{$nNameIndex}{FA} = $normal_FA;
		$att->{$tNameIndex}{FA} = $tumor_FA;
	}
	
	### Check for GT
	unless (exists $att->{$nNameIndex}{GT}) {
		# stupid strelka doesn't put in genotype information
		# since GT is usually first in attributes, we'll put it there
		# normal genotype use the info NT key
		if ($att->{INFO}{NT} =~ /ref/i) {
			$att->{$nNameIndex}{GT} = '0/0';
		}
		elsif ($att->{INFO}{NT} =~ /het/i) {
			# info value SGT=ref->hom
			$att->{$nNameIndex}{GT} = '0/1';
		}
		elsif ($att->{INFO}{NT} =~ /hom$/i) {
			# info value SGT=ref->hom
			$att->{$nNameIndex}{GT} = '1/1';
		}
		else {
			warn sprintf "no Normal genotype key NT for data line %s:%d\n", 
				$line->seq_id, $line->start;
			$att->{$nNameIndex}{GT} = './.';
		}
		
		# tumor genotype uses the info key SGT
		if ($att->{INFO}{SGT} =~ /het$/i) {
			# info value SGT=ref->het;
			$att->{$tNameIndex}{GT} = '0/1';
		}
		elsif ($att->{INFO}{SGT} =~ /hom$/i) {
			# info value SGT=ref->hom
			$att->{$tNameIndex}{GT} = '1/1';
		}
		else {
			warn sprintf "no Normal genotype key SGT for data line %s:%d\n", 
				$line->seq_id, $line->start;
			$att->{$tNameIndex}{GT} = './.';
		}
	}
	
	### Check for depth
	if (not exists $att->{$nNameIndex}{AD} and not exists $att->{$nNameIndex}{DP4}) {
		# this should only be for the Strelka stuff
		$tier ||= 1; # in case it wasn't defined for some reason
		
		# reference count
		my $normalRefCount = $tier == 1 ? $att->{$nNameIndex}{DP} : $att->{$nNameIndex}{DP2};
		my $tumorRefCount  = $tier == 1 ? $att->{$tNameIndex}{DP} : $att->{$tNameIndex}{DP2};
		
		# alternate count
		my ($normalAltCount, $tumorAltCount);
			# we'll need to split the appropriate FORMAT attribute into tier1,tier2
			# then take the appropriate index for the tier
		if ($indel) {
			$normalAltCount = (split ',', $att->{$nNameIndex}{TIR}, 2)[ $tier - 1 ];
			$tumorAltCount  = (split ',', $att->{$tNameIndex}{TIR}, 2)[ $tier - 1 ];
		} 
		else {
			$normalAltCount = (split ',', $att->{$nNameIndex}{TAR}, 2)[ $tier - 1 ];
			$tumorAltCount  = (split ',', $att->{$tNameIndex}{TAR}, 2)[ $tier - 1 ];
		}
		
		# put back
		$att->{$nNameIndex}{AD} = "$normalRefCount,$normalAltCount";
		$att->{$tNameIndex}{AD} = "$tumorRefCount,$tumorAltCount";
	}
	elsif (not exists $att->{$nNameIndex}{AD} and exists $att->{$nNameIndex}{DP4}) {
		# this should be for SomaticSniper 
		# this is a temporary fix, probably should not be kept permanently
		# convert AD ref,alt from DP4 
		my @nC = split ',', $att->{$nNameIndex}{DP4}, 4;
		my @tC = split ',', $att->{$tNameIndex}{DP4}, 4;
		$att->{$nNameIndex}{AD} = sprintf("%d,%d", $nC[0] + $nC[1], $nC[2] + $nC[3]);
		$att->{$tNameIndex}{AD} = sprintf("%d,%d", $tC[0] + $tC[1], $tC[2] + $tC[3]);
	}
	elsif (exists $att->{$nNameIndex}{AD} and 
		exists $att->{$nNameIndex}{RD} and 
		exists $att->{$nNameIndex}{DP4}
	) {
		# VarScan is using separate counts for everything
		$att->{$nNameIndex}{AD} = sprintf("%d,%d", $att->{$nNameIndex}{RD}, 
			$att->{$nNameIndex}{AD});
		$att->{$tNameIndex}{AD} = sprintf("%d,%d", $att->{$tNameIndex}{RD}, 
			$att->{$tNameIndex}{AD});
		$change_AD_comment++;
	}
	
	### update the line
	if ($minimize) {
		# we need to delete extraneous attributes
		
		# sample columns
		foreach (keys %{$att->{$nNameIndex}}) {
			next if ($_ eq 'GT' or $_ eq 'AD' or $_ eq 'FA');
			delete $att->{$nNameIndex}{$_};
			delete $att->{$tNameIndex}{$_};
		}
		
		# INFO column
		foreach (keys %{$att->{INFO}} ) {
			# basically remove everything in the INFO field
			# the idea is to minimize, right????
			delete $att->{INFO}{$_};
		}
	}
	$line->rewrite_vcf_attributes or warn "failed to rewrite vcf attributes!\n";
}



# update the VCF metadata comment lines as necessary
my $vcf_head = $Data->vcf_headers;
if ($minimize) {
	# delete unnecessary FORMAT keys
	foreach my $key (keys %{ $vcf_head->{FORMAT} }) {
		next if ($key eq 'GT' or $key eq 'AD' or $key eq 'FA');
		delete $vcf_head->{FORMAT}{$key};
	}

	# delete unnecessary INFO keys
	foreach (keys %{$vcf_head->{INFO}}) {
		delete $vcf_head->{INFO}{$_};
	}
}
	
unless (exists $vcf_head->{FORMAT}{GT}) {
	# add FA metadata line
	$vcf_head->{FORMAT}{GT} = q(ID=GT,Number=1,Type=String,Description="Genotype");
}
unless (exists $vcf_head->{FORMAT}{FA}) {
	# add FA metadata line
	$vcf_head->{FORMAT}{FA} = q(<ID=FA,Number=A,Type=Float,Description="Allele fraction of the alternate allele with regard to reference">);
}
if ($change_AD_comment or not exists $vcf_head->{FORMAT}{AD}) {
	# add FA metadata line
	$vcf_head->{FORMAT}{AD} = q(ID=AD,Number=.,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed");
}
	
# rewrite the VCF metadata headers
$Data->rewrite_vcf_headers;



# write the file
# do not gzip, user should use bgzip and index themselves
my $success = $Data->write_file(filename => $outfile, gz => 0);
print "Wrote updated $success\n";

__END__

=head2 MuTect FORMAT keys

	##FORMAT=<ID=AD,Number=.,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">
	##FORMAT=<ID=BQ,Number=A,Type=Float,Description="Average base quality for reads supporting alleles">
	##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth (reads with MQ=255 or with bad mates are filtered)">
	##FORMAT=<ID=FA,Number=A,Type=Float,Description="Allele fraction of the alternate allele with regard to reference">
	##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
	##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
	##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification">
	##FORMAT=<ID=SS,Number=1,Type=Integer,Description="Variant status relative to non-adjacent Normal,0=wildtype,1=germline,2=somatic,3=LOH,4=post-transcriptional modification,5=unknown">

=head2 SomaticSniper FORMAT keys

	##FORMAT=<ID=AMQ,Number=.,Type=Integer,Description="Average mapping quality for each allele present in the genotype">
	##FORMAT=<ID=BCOUNT,Number=4,Type=Integer,Description="Occurrence count for each base at this site (A,C,G,T)">
	##FORMAT=<ID=BQ,Number=.,Type=Integer,Description="Average base quality">
	##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Total read depth">
	##FORMAT=<ID=DP4,Number=4,Type=Integer,Description="# high-quality ref-forward bases, ref-reverse, alt-forward and alt-reverse bases">
	##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">
	##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
	##FORMAT=<ID=IGT,Number=1,Type=String,Description="Genotype when called independently (only filled if called in joint prior mode)">
	##FORMAT=<ID=JGQ,Number=1,Type=Integer,Description="Joint genotype quality (only filled if called in join prior mode)">
	##FORMAT=<ID=MQ,Number=1,Type=Integer,Description="Average mapping quality across all reads">
	##FORMAT=<ID=SS,Number=1,Type=Integer,Description="Variant status relative to non-adjacent Normal, 0=wildtype,1=germline,2=somatic,3=LOH,4=unknown">
	##FORMAT=<ID=SSC,Number=1,Type=Integer,Description="Somatic Score">
	##FORMAT=<ID=VAQ,Number=1,Type=Integer,Description="Variant allele quality">

=head2 Strelka Indel FORMAT keys

	##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth for tier1">
	##FORMAT=<ID=DP2,Number=1,Type=Integer,Description="Read depth for tier2">
	##FORMAT=<ID=DP50,Number=1,Type=Float,Description="Average tier1 read depth within 50 bases">
	##FORMAT=<ID=FDP50,Number=1,Type=Float,Description="Average tier1 number of basecalls filtered from original read depth within 50 bases">
	##FORMAT=<ID=SUBDP50,Number=1,Type=Float,Description="Average number of reads below tier1 mapping quality threshold aligned across sites within 50 bases">
	##FORMAT=<ID=TAR,Number=2,Type=Integer,Description="Reads strongly supporting alternate allele for tiers 1,2">
	##FORMAT=<ID=TIR,Number=2,Type=Integer,Description="Reads strongly supporting indel allele for tiers 1,2">
	##FORMAT=<ID=TOR,Number=2,Type=Integer,Description="Other reads (weak support or insufficient indel breakpoint overlap) for tiers 1,2">

=head2 VarScan FORMAT keys

	##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
	##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
	##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">
	##FORMAT=<ID=RD,Number=1,Type=Integer,Description="Depth of reference-supporting bases (reads1)">
	##FORMAT=<ID=AD,Number=1,Type=Integer,Description="Depth of variant-supporting bases (reads2)">
	##FORMAT=<ID=FREQ,Number=1,Type=String,Description="Variant allele frequency">
	##FORMAT=<ID=DP4,Number=1,Type=String,Description="Strand read counts: ref/fwd, ref/rev, var/fwd, var/rev">


-head2 Scalpel FORMAT keys
	
	##FORMAT=<ID=AD,Number=.,Type=Integer,Description="k-mer depth supporting reference/indel at the site">
	##FORMAT=<ID=DP,Number=1,Type=Integer,Description="k-mer Depth">
	##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">


=head2 ideal minimal
	GT
	AD (DP4)
	FA


