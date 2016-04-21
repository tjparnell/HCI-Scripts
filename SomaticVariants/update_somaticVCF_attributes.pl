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
use Bio::ToolBox::Data 1.34;
my $VERSION = 1.2;

unless (@ARGV) {
	print <<END;

A script to fix and standardize sample attributes in somatic VCF files. 
  - adds missing alternate fraction (FA or FREQ) attribute (Strelka and SomaticSniper)
  - adds missing genotype GT attribute (Strelka)
  - adds missing allele depth (AD or DP4) attribute (Strelka)
  - optionally discard excess sample attributes, leaving only GT:AD:FA
      otherwise merging VCFs from different callers with different sample attributes 
      into a single VCF often leads to malformed VCF sample records

Usage: $0 -i <input.vcf> 
Options: 
  -o <output.vcf>       default is to overwrite input!!!!
  -t <NAME>             Tumor sample name, default TUMOR
  -n <NAME>             Normal sample name, default NORMAL
  -m                    flag to minimize to just three sample attributes: GT:AD:FA
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
my @noINFO = qw(IC IHP NT OVERLAP QSI QSI_NT RC RU SGT SVTYPE TQSI TQSI_NT MQ0 SOMATIC VT DB .);

# input file
my $Data = Bio::ToolBox::Data->new(in => $file) or 
	die "unable to open $file!\n";
die "file is not VCF!!!\n" unless $Data->vcf;

# output file
$outfile ||= $file;
$outfile =~ s/\.gz$//i;
$outfile = $outfile . '.vcf' unless $outfile =~ /\.vcf$/i;

# sample columns
my $norm_i = $Data->find_column($nName);
my $tumr_i = $Data->find_column($tName);
unless ($norm_i) {
	die "can not find index for $nName!\n";
}
unless ($tumr_i) {
	die "can not find index for $tName!\n";
}

# process
my $iterator = $Data->row_stream;
while (my $line = $iterator->next_row) {
	
	my $att = $line->vcf_attributes;
	
	# Strelka data tier and variant type, used in both FA and depth calculations
	my ($tier, $indel);
	
	### check for FA or FREQ
	unless (exists $att->{$nName}{FA} or exists $att->{$nName}{FREQ}) {
		# No frequency of alternate
		
		# calculate FREQ attribute
		my ($normal_FA, $tumor_FA);
		if (exists $att->{$nName}{DP4}) {
			# DP4 = ref-forward bases, ref-reverse, alt-forward and alt-reverse bases
			my ($rf, $rr, $af, $ar) = split ',', $att->{$nName}{DP4}, 4;
			$normal_FA = sprintf "%0.2f", ($af + $ar) / ($rf + $rr + $af + $ar);
			($rf, $rr, $af, $ar) = split ',', $att->{$tName}{DP4}, 4;
			$tumor_FA = sprintf "%0.2f", ($af + $ar) / ($rf + $rr + $af + $ar);
		} 
		
		elsif (exists $att->{$nName}{TIR} or exists $att->{$nName}{TAR}) {
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
								$indel ? $att->{$nName}{TIR} : $att->{$nName}{TAR}, 
								2;
			my ($tt1, $tt2) = split ',', # tumor
								$indel ? $att->{$tName}{TIR} : $att->{$tName}{TAR}, 
								2;
			
			# calculate FREQ based on which data tier was used
			# using the total depth DP as the denominator
			if ($tier == 1) {
				$normal_FA = $att->{$nName}{DP} ? 
					sprintf "%0.2f", $nt1 / $att->{$nName}{DP} : 0;
				$tumor_FA  = $att->{$tName}{DP} ? 
					sprintf "%0.2f", $tt1 / $att->{$tName}{DP} : 0;
			} elsif ($tier == 2) {
				$normal_FA = $att->{$nName}{DP2} ? 
					sprintf "%0.2f", $nt1 / $att->{$nName}{DP2} : 0;
				$tumor_FA  = $att->{$tName}{DP2} ? 
					sprintf "%0.2f", $tt1 / $att->{$tName}{DP2} : 0;
			} else {
				warn sprintf "unknown tier '$tier' for data line %s:%d\n",
					$line->seq_id, $line->start;
				$normal_FA = 0;
				$tumor_FA  = 0;
			}
		} 
		
		elsif (exists $att->{$nName}{AD}) {
			# AD = ref bases, alt bases
			my ($ref, $alt) = split ',', $att->{$nName}{AD}, 2;
			$normal_FA = sprintf "%0.2f", $alt / ($ref + $alt);
			($ref, $alt) = split ',', $att->{$tName}{AD}, 2;
			$tumor_FA = sprintf "%0.2f", $alt / ($ref + $alt);
		} 
		
		else {
			warn "missing read depth identifiers for data line %s:%d\n",
				$line->seq_id, $line->start;
			$normal_FA = 0;
			$tumor_FA  = 0;
		}
		
		# push the new FA attribute
		$att->{$nName}{FA} = $normal_FA;
		$att->{$tName}{FA} = $tumor_FA;
	}
	
	### Check for GT
	unless (exists $att->{$nName}{GT}) {
		# stupid strelka doesn't put in genotype information
		# since GT is usually first in attributes, we'll put it there
		# normal genotype use the info NT key
		if ($att->{INFO}{NT} =~ /ref/i) {
			$att->{$nName}{GT} = '0/0';
		}
		elsif ($att->{INFO}{NT} =~ /het/i) {
			# info value SGT=ref->hom
			$att->{$nName}{GT} = '0/1';
		}
		elsif ($att->{INFO}{NT} =~ /hom$/i) {
			# info value SGT=ref->hom
			$att->{$nName}{GT} = '1/1';
		}
		else {
			warn sprintf "no Normal genotype key NT for data line %s:%d\n", 
				$line->seq_id, $line->start;
			$att->{$nName}{GT} = './.';
		}
		
		# tumor genotype uses the info key SGT
		if ($att->{INFO}{SGT} =~ /het$/i) {
			# info value SGT=ref->het;
			$att->{$tName}{GT} = '0/1';
		}
		elsif ($att->{INFO}{SGT} =~ /hom$/i) {
			# info value SGT=ref->hom
			$att->{$tName}{GT} = '1/1';
		}
		else {
			warn sprintf "no Normal genotype key SGT for data line %s:%d\n", 
				$line->seq_id, $line->start;
			$att->{$tName}{GT} = './.';
		}
	}
	
	### Check for depth
	if (not exists $att->{$nName}{AD} and not exists $att->{$nName}{DP4}) {
		# this should only be for the Strelka stuff
		$tier ||= 1; # in case it wasn't defined for some reason
		
		# reference count
		my $normalRefCount = $tier == 1 ? $att->{$nName}{DP} : $att->{$nName}{DP2};
		my $tumorRefCount  = $tier == 1 ? $att->{$tName}{DP} : $att->{$tName}{DP2};
		
		# alternate count
		my ($normalAltCount, $tumorAltCount);
			# we'll need to split the appropriate FORMAT attribute into tier1,tier2
			# then take the appropriate index for the tier
		if ($indel) {
			$normalAltCount = (split ',', $att->{$nName}{TIR}, 2)[ $tier - 1 ];
			$tumorAltCount  = (split ',', $att->{$tName}{TIR}, 2)[ $tier - 1 ];
		} else {
			$normalAltCount = (split ',', $att->{$nName}{TAR}, 2)[ $tier - 1 ];
			$tumorAltCount  = (split ',', $att->{$tName}{TAR}, 2)[ $tier - 1 ];
		}
		
		# put back
		$att->{$nName}{AD} = "$normalRefCount,$normalAltCount";
		$att->{$tName}{AD} = "$tumorRefCount,$tumorAltCount";
	}
	elsif (not exists $att->{$nName}{AD} and exists $att->{$nName}{DP4}) {
		# this should be for SomaticSniper 
		# this is a temporary fix, probably should not be kept permanently
		# convert AD ref,alt from DP4 
		my @nC = split ',', $att->{$nName}{DP4}, 4;
		my @tC = split ',', $att->{$tName}{DP4}, 4;
		my $normalCount = sprintf("%d,%d", $nC[0] + $nC[1], $nC[2] + $nC[3]);
		my $tumorCount = sprintf("%d,%d", $tC[0] + $tC[1], $tC[2] + $tC[3]);
		
		$att->{$nName}{AD} = $normalCount;
		$att->{$tName}{AD} = $tumorCount;
	}
	
	
	
	# update the line
	if ($minimize) {
		
		# format column
		$line->value(8, join(':', qw(GT AD FA)) );
		
		# sample columns
		$line->value($norm_i, join(':', map {$att->{$nName}{$_}} qw(GT AD FA)) );
		$line->value($tumr_i, join(':', map {$att->{$tName}{$_}} qw(GT AD FA)) );
		
		# update INFO column
		foreach my $id (@noINFO) {
			# these are the extraneous fields found in the somatic callers I've been using
			# keep anything that's not in this list
			if (exists $att->{INFO}{$id}) {
				delete $att->{INFO}{$id};
			}
		}
		
		# regenerate INFO field with remainder
		my $info = join(';', 
			map { join('=', $_, $att->{INFO}{$_}) } 
			sort {$a cmp $b} 
			keys %{$att->{INFO}}
		);
		$info ||= '.'; # sometimes we have nothing left
		$line->value(7, $info);
	}
	else {
		# we can skip the INFO column since we don't modify that
		
		# determine the format order, alphabetize since we can't maintain original order
		my @formatKeys = qw(GT); # genotype must come first
		foreach (sort {$a cmp $b} keys %{ $att->{$nName} }) {
			push @formatKeys, $_ unless $_ eq 'GT';
		}
		$line->value(8, join(":", @formatKeys));
		$line->value($norm_i, join(":", map { $att->{$nName}{$_} } @formatKeys ) );
		$line->value($tumr_i, join(":", map { $att->{$tName}{$_} } @formatKeys ) );
	}
}

# update the VCF metadata comment lines as necessary
if ($minimize) {
	my @comments = $Data->comments; # returns an array
	my @to_delete;
	
	# format metadata lines
	foreach (my $i = 0; $i < scalar @comments; $i++) {
		if ($comments[$i] =~ /^##FORMAT/) {
			push @to_delete, $i unless $comments[$i] =~ /ID=(?:GT|AD|FA)/;
		}
	}
	
	# info metadata lines
	foreach my $id (@noINFO) {
		# these are the extraneous fields found in the somatic callers I've been using
		# keep anything that's not in this list
		foreach (my $i = 0; $i < scalar @comments; $i++) {
			if ($comments[$i] =~ /^##INFO/) {
				push @to_delete, $i if $comments[$i] =~ /ID=$id/;
			}
		}
	}
	
	# update VCF metadata comment lines
	foreach (sort {$b <=> $a} @to_delete) {
		# must reverse sort to avoid messing up other lines
		$Data->delete_comment($_);
	}
}



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


=head2 ideal minimal
	GT
	AD (DP4)
	FA


