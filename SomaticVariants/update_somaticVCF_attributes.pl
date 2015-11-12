#!/usr/bin/env perl

use strict;
use Getopt::Long;
# use Data::Dumper;
use Bio::ToolBox::Data;
	# the new Bio::ToolBox version has a VCF attributes method call, but 
	# this is older script, so we're getting the attributes manually here

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
  -n <NAME>             Normal sample name, default NORMAL
  -t <NAME>             Tumor sample name, default TUMOR
  -m                    flag to minimize to just three sample attributes: GT:AD:FA
END
	exit;
}

my $file;
my $outfile;
my $normName = 'NORMAL';
my $tumorName = 'TUMOR';
my $minimize;

GetOptions( 
	'i=s'       => \$file, # input vcf
	'o=s'       => \$outfile, # output vcf
	'n=s'       => \$normName, # normal name
	't=s'       => \$tumorName, # tumor name
	'm!'        => \$minimize, # throw away extra attributes
) or die "bad options!\n";


my $Data = Bio::ToolBox::Data->new(in => $file) or 
	die "unable to open $file!\n";
$outfile ||= $file;
$outfile =~ s/\.gz$//i;
$outfile = $outfile . '.vcf' unless $outfile =~ /\.vcf$/i;

my $info_i = $Data->find_column('INFO');
my $form_i = $Data->find_column('FORMAT');
my $norm_i = $Data->find_column($normName);
my $tumr_i = $Data->find_column($tumorName);
unless ($norm_i) {
	die "can not find index for $normName!\n";
}
unless ($tumr_i) {
	die "can not find index for $tumorName!\n";
}

my $iterator = $Data->row_stream;
while (my $line = $iterator->next_row) {
	
	# split values
	my @formatKeys = split /:/, $line->value($form_i);
	my @normalVals = split /:/, $line->value($norm_i);
	my @tumorVals  = split /:/, $line->value($tumr_i);
	my $number = scalar @formatKeys;
	warn "format and sample attributes have unequal number of entries for data line " . 
		$line->seq_id . ":" . $line->start . "\n" if 
		(scalar @normalVals != $number or scalar @tumorVals != $number);
	
	# generate attribute hash for each sample
	# this could be updated with new version of Bio::ToolBox to do this automatically
	my %normals = map { $formatKeys[$_] => $normalVals[$_] } (0 .. $#formatKeys);
	my %tumors  = map { $formatKeys[$_] => $tumorVals[$_] } (0 .. $#formatKeys);
	my %info = 	map { $_->[0] => defined $_->[1] ? $_->[1] : 1 } 
				map { [split(/=/, $_)] } 
				split(/;/, $line->value($info_i));
# 	printf "Data for variant %s:%s\nINFO %s\nFORMAT %s\nNORMAL %s\nTUMOR %s\nNORMALHASH %s\nTUMORHASH %s\n", 
# 		$line->seq_id, $line->start, Dumper(\%info), Dumper(\@formatKeys), 
# 		Dumper(\@normalVals), Dumper(\@tumorVals), Dumper(\%normals), Dumper(\%tumors);
	
	# Strelka data tier and variant type, used in both FA and depth calculations
	my ($tier, $indel);
	
	### check for FA or FREQ
	unless (exists $normals{FA} or exists $normals{FREQ}) {
		# No frequency of alternate
		
		# calculate FREQ attribute
		my ($normal_FA, $tumor_FA);
		if (exists $normals{DP4}) {
			# DP4 = ref-forward bases, ref-reverse, alt-forward and alt-reverse bases
			my ($rf, $rr, $af, $ar) = split /,/, $normals{DP4};
			$normal_FA = sprintf "%0.2f", ($af + $ar) / ($rf + $rr + $af + $ar);
			($rf, $rr, $af, $ar) = split /,/, $tumors{DP4};
			$tumor_FA = sprintf "%0.2f", ($af + $ar) / ($rf + $rr + $af + $ar);
		} 
		
		elsif (exists $normals{TIR} or exists $normals{TAR}) {
			# Strelka attributes
			
			# determine data tier and type of indel
			if (exists $info{TQSI}) {
				# strelka indel
				$tier = $info{TQSI};
				$indel = 1;
			}
			elsif (exists $info{TQSS}) {
				# strelka somatic
				$tier = $info{TQSS};
				$indel = 0;
			}
			unless ($tier) {
				# just go with the first data tier 
				warn "no tier data for data line " . $line->seq_id, ":" , $line->start, "\n";
				$tier = 1;
			}
			
			# data tier read counts appropriate for somatic variant type
			my ($nt1, $nt2) = split /,/, $indel ? $normals{TIR} : $normals{TAR};
			my ($tt1, $tt2) = split /,/, $indel ? $tumors{TIR}  : $tumors{TAR};
			
			# calculate FREQ based on which data tier was used
			if ($tier == 1) {
				$normal_FA = $normals{DP} ? sprintf "%0.2f", $nt1 / $normals{DP} : 0;
				$tumor_FA  = $tumors{DP} ? sprintf "%0.2f", $tt1 / $tumors{DP} : 0;
			} elsif ($tier == 2) {
				$normal_FA = $normals{DP2} ? sprintf "%0.2f", $nt2 / $normals{DP2} : 0;
				$tumor_FA  = $tumors{DP2} ? sprintf "%0.2f", $tt2 / $tumors{DP2} : 0;
			} else {
				warn "unknown tier '$tier' for data line " . $line->seq_id, ":" , $line->start, "\n";
				$normal_FA = 0;
				$tumor_FA  = 0;
			}
		} else {
			warn "missing read depth identifiers for data line " . $line->seq_id, ":" , $line->start, "\n";
			$normal_FA = 0;
			$tumor_FA  = 0;
		}
		
		# push the new FA attribute
		$normals{FA} =  $normal_FA;
		$tumors{FA} = $tumor_FA;
		push @formatKeys, 'FA';
		push @normalVals, $normal_FA;
		push @tumorVals, $tumor_FA;
	}
	
	### Check for GT
	unless (exists $normals{GT}) {
		# stupid strelka doesn't put in genotype information
		# since GT is usually first in attributes, we'll put it there
		unshift @formatKeys, 'GT';
		# normal genotype use the info NT key
		if ($info{NT} =~ /ref/i) {
			unshift @normalVals, '0/0';
			$normals{GT} = '0/0';
		}
		elsif ($info{NT} =~ /het/i) {
			# info value SGT=ref->hom
			unshift @normalVals, '0/1';
			$normals{GT} = '0/1';
		}
		elsif ($info{NT} =~ /hom$/i) {
			# info value SGT=ref->hom
			unshift @normalVals, '1/1';
			$normals{GT} = '1/1';
		}
		else {
			warn "no Normal genotype key NT for data line " . $line->seq_id, ":" , $line->start, "\n";
			unshift @normalVals, './.';
			$normals{GT} = './.';
		}
		
		# tumor genotype uses the info key SGT
		if ($info{SGT} =~ /het$/i) {
			# info value SGT=ref->het;
			unshift @tumorVals, '0/1';
			$tumors{GT} = '0/1';
		}
		elsif ($info{SGT} =~ /hom$/i) {
			# info value SGT=ref->hom
			unshift @tumorVals, '1/1';
			$tumors{GT} = '1/1';
		}
		else {
			warn "no Tumor genotype key SGT for data line " . $line->seq_id, ":" , $line->start, "\n";
			unshift @tumorVals, './.';
			$tumors{GT} = './.';
		}
	}
	
	### Check for depth
	if (not exists $normals{AD} and not exists $normals{DP4}) {
		# this should only be for the Strelka stuff
		$tier ||= 1; # in case it wasn't defined for some reason
		
		# reference count
		my $normalRefCount = $tier == 1 ? $normals{DP} : $normals{DP2};
		my $tumorRefCount  = $tier == 1 ? $tumors{DP}  : $tumors{DP2};
		
		# alternate count
		my ($normalAltCount, $tumorAltCount);
			# we'll need to split the appropriate FORMAT attribute into tier1,tier2
			# then take the appropriate index for the tier
		if ($indel) {
			$normalAltCount = (split ',', $normals{TIR})[ $tier - 1 ];
			$tumorAltCount  = (split ',', $tumors{TIR})[ $tier - 1 ];
		} else {
			$normalAltCount = (split ',', $normals{TAR})[ $tier - 1 ];
			$tumorAltCount  = (split ',', $tumors{TAR})[ $tier - 1 ];
		}
		
		# put back
		$normals{AD} = "$normalRefCount,$normalAltCount";
		$tumors{AD} = "$tumorRefCount,$tumorAltCount";
		push @formatKeys, 'AD';
		push @normalVals, "$normalRefCount,$normalAltCount";
		push @tumorVals, "$tumorRefCount,$tumorAltCount";
	}
	elsif (not exists $normals{AD} and exists $normals{DP4}) {
		# this should be for SomaticSniper 
		# this is a temporary fix, probably should not be kept permanently
		# convert AD ref,alt from DP4 
		my @nC = split ",", $normals{DP4};
		my @tC = split ",", $tumors{DP4};
		my $normalCount = sprintf("%s,%s", $nC[0] + $nC[1], $nC[2] + $nC[3]);
		my $tumorCount = sprintf("%s,%s", $tC[0] + $tC[1], $tC[2] + $tC[3]);
		
		$normals{AD} = $normalCount;
		$tumors{AD} = $tumorCount;
		push @formatKeys, 'AD';
		push @normalVals, $normalCount;
		push @tumorVals, $tumorCount;
	}
	
	### Check for depth again
	
	
	
	# update the line
	if ($minimize) {
		$line->value($form_i, join(':', qw(GT AD FA)) );
		$line->value($norm_i, join(':', map {$normals{$_}} qw(GT AD FA)) );
		$line->value($tumr_i, join(':', map {$tumors{$_}} qw(GT AD FA)) );
	}
	else {
		$line->value($form_i, join(":", @formatKeys));
		$line->value($norm_i, join(":", @normalVals));
		$line->value($tumr_i, join(":", @tumorVals));
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


