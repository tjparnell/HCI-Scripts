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
use Bio::ToolBox::Data '1.41';
my $VERSION = 1.2;

unless (scalar @ARGV) {
	print <<END;

A script to convert other style hg38 chromosome identifiers to Ensembl identifiers.
This uses built-in synonyms for UCSC, NCBI, and others.

Usage:
 $0 <infile> <outfile>

It will handle BED, GFF, GTF, refFlat, genePred, and VCF files.

It will also handle SAM header files, which can be used in 
conjunction with samtools reheader to rename Bam files.

It will report any unmatched chromosomes that couldn't be converted. 
These are skipped in the output file.

END
	exit;
}

my $infile = shift @ARGV;
my $outfile = shift @ARGV || undef;
my $lookup = make_lookup_table();

if ($infile =~ /\.fa(?:sta)?(?:\.gz)?$/i) {
	# input is a fasta file
	die "can't do fasta files!!!"
}
elsif ($infile =~ /\.sam$/) {
	# a sam file
	process_sam();
}
else {
	# otherwise assume some sort of gene table
	process_table();
}

sub make_lookup_table {
	my %lookup;
	while (my $line = <DATA>) {
		chomp $line;
		my ($ucsc, $ens) = split('\t', $line);
		$lookup{$ucsc} = $ens;
	}
	return \%lookup;
}

sub process_table {
	my $Stream = Bio::ToolBox::Data->new(
		stream => 1,
		in     => $infile,
	) or die " unable to open input table $infile! $!\n";
	# we check the chromosome column below
	
	# make output
	unless ($outfile) {
		$outfile = $Stream->path . $Stream->basename . '.ensembl' . $Stream->extension;
	}
	my $Out = $Stream->duplicate($outfile) or 
		die " unable to open output stream for file $outfile! $!\n";
	
	# deal with metadata
	my @comments = $Stream->comments;
	if (@comments) {
		for (my $i =$#comments; $i >= 0; $i--) {
			# delete the existing comments, these are indexed so go in reverse
			# order, we'll add back fixed ones
			$Out->delete_comment($i); 
		}
		foreach my $c (@comments) {
			if ($c =~ /^##sequence\-region\s+([\w\.]+)\s/) {
				# gff3 sequence pragmas
				my $chr = $1;
				if (exists $lookup->{$chr}) {
					my $alt = $lookup->{$chr};
					$c =~ s/$chr/$alt/;
				}
			}
			elsif ($c =~ /^##contig=<ID=([\w\.]+)/) {
				# vcf sequence identifiers
				my $chr = $1;
				if (exists $lookup->{$chr}) {
					my $alt = $lookup->{$chr};
					$c =~ s/$chr/$alt/;
				}
			}
			$Out->add_comment($c);
		}
		
	}
	
	# data replacements
	my $seq_i = $Stream->chromo_column;
	die "can't find chromosome column!\n" unless defined $seq_i;
	my %notfound;
	while (my $row = $Stream->next_row) {
		my $chr = $row->value($seq_i);
		if (exists $lookup->{$chr}) {
			$row->value($seq_i, $lookup->{$chr});
			$Out->write_row($row);
		}
		else {
			$notfound{$chr}++;
		}
	}
	$Out->close_fh;
	$Stream->close_fh;
	if (%notfound) { 
		printf "could not convert the following chromosomes:\n%s\n", 
			join("\n", map {sprintf(" %-25s%d lines", $_, $notfound{$_})} keys %notfound);
	}
	printf "wrote %s\n", $Out->filename;
}


sub process_sam {
	# open input
	my $in_fh = Bio::ToolBox::Data->open_to_read_fh($infile) 
		or die "unable to open '$infile'! $!\n";
	
	
	# make output
	unless ($outfile) {
		$outfile = $infile;
		$outfile =~ s/sam$/ensembl.sam/i;
	}
	my $out_fh = Bio::ToolBox::Data->open_to_write_fh($outfile) or 
		die " unable to open output file handle for $outfile! $!\n";
	
	# process
	my %notfound;
	while (my $line = $in_fh->getline) {
		# check header
		if (substr($line,0,3) eq q(@SQ)) {
			# a sequence header line
			my @bits = split("\t", $line);
			if ($bits[1] =~ /^SN:([\w\.]+)$/) {
				if (exists $lookup->{$1}) {
					$bits[1] = sprintf("SN:%s", $lookup->{$1});
				}
				else {
					$notfound{$1}++;
				}
			}
			else {
				# not a classic order with SN first!????
				warn "non-standard \@SQ line!\n";
			}
			# replace
			$line = join("\t", @bits);
		}
		elsif (substr($line,0,1) eq q(@)) {
			# some other header line
		}
		else {
			# alignment line!?
			my @bits = split("\t", $line);
			if ($bits[2] ne '*') {
				if (exists $lookup->{$bits[2]}) {
					$bits[2] = $lookup->{ $bits[2] };
				}
				else {
					$notfound{$bits[2]}++;
				}
			}
			if ($bits[6] ne '*' or $bits[6] ne '=') {
				if (exists $lookup->{$bits[6]}) {
					$bits[6] = $lookup->{ $bits[6] };
				}
				else {
					$notfound{$bits[6]}++;
				}
			}
			$line = join("\t", @bits);
		}
		
		# print the line
		$out_fh->print($line);
	}
	
	# Finished
	$out_fh->close;
	$in_fh->close;
	if (%notfound) { 
		printf "could not convert the following chromosomes:\n%s\n", 
			join("\n", map {sprintf(" %-25s%d lines", $_, $notfound{$_})} keys %notfound);
	}
	printf "wrote %s\n", $outfile;
}


# Manually curated from UCSC hg38 and Ensembl GRCh38.dna.primary_assembly
# does not include alt chromsomes at this time
__DATA__
chr1	1
chr2	2
chr3	3
chr4	4
chr5	5
chr6	6
chr7	7
chr8	8
chr9	9
chr11	11
chr10	10
chr12	12
chr13	13
chr14	14
chr15	15
chr16	16
chr17	17
chr18	18
chr20	20
chr19	19
chr22	22
chr21	21
chrX	X
chrY	Y
chrM	MT
chr16_KI270728v1_random	KI270728.1
chr15_KI270727v1_random	KI270727.1
chrUn_KI270442v1	KI270442.1
chr17_KI270729v1_random	KI270729.1
chr14_GL000225v1_random	GL000225.1
chrUn_KI270743v1	KI270743.1
chr4_GL000008v2_random	GL000008.2
chr14_GL000009v2_random	GL000009.2
chrUn_KI270747v1	KI270747.1
chr14_KI270722v1_random	KI270722.1
chr14_GL000194v1_random	GL000194.1
chrUn_KI270742v1	KI270742.1
chr17_GL000205v2_random	GL000205.2
chrUn_GL000195v1	GL000195.1
chr22_KI270736v1_random	KI270736.1
chr22_KI270733v1_random	KI270733.1
chrUn_GL000224v1	GL000224.1
chrUn_GL000219v1	GL000219.1
chr9_KI270719v1_random	KI270719.1
chrUn_GL000216v2	GL000216.2
chr1_KI270712v1_random	KI270712.1
chr1_KI270706v1_random	KI270706.1
chr14_KI270725v1_random	KI270725.1
chrUn_KI270744v1	KI270744.1
chr22_KI270734v1_random	KI270734.1
chrUn_GL000213v1	GL000213.1
chrUn_GL000220v1	GL000220.1
chr2_KI270715v1_random	KI270715.1
chrUn_GL000218v1	GL000218.1
chrUn_KI270749v1	KI270749.1
chrUn_KI270741v1	KI270741.1
chr3_GL000221v1_random	GL000221.1
chr2_KI270716v1_random	KI270716.1
chr22_KI270731v1_random	KI270731.1
chrUn_KI270751v1	KI270751.1
chrUn_KI270750v1	KI270750.1
chrUn_KI270519v1	KI270519.1
chrUn_GL000214v1	GL000214.1
chr1_KI270708v1_random	KI270708.1
chr17_KI270730v1_random	KI270730.1
chrUn_KI270438v1	KI270438.1
chr22_KI270737v1_random	KI270737.1
chr11_KI270721v1_random	KI270721.1
chr22_KI270738v1_random	KI270738.1
chrUn_KI270748v1	KI270748.1
chrUn_KI270435v1	KI270435.1
chr5_GL000208v1_random	GL000208.1
chrUn_KI270538v1	KI270538.1
chrUn_KI270756v1	KI270756.1
chr22_KI270739v1_random	KI270739.1
chrUn_KI270757v1	KI270757.1
chr1_KI270709v1_random	KI270709.1
chrUn_KI270746v1	KI270746.1
chrUn_KI270753v1	KI270753.1
chrUn_KI270589v1	KI270589.1
chr14_KI270726v1_random	KI270726.1
chr22_KI270735v1_random	KI270735.1
chr1_KI270711v1_random	KI270711.1
chrUn_KI270745v1	KI270745.1
chr1_KI270714v1_random	KI270714.1
chr22_KI270732v1_random	KI270732.1
chr1_KI270713v1_random	KI270713.1
chrUn_KI270754v1	KI270754.1
chr1_KI270710v1_random	KI270710.1
chr9_KI270717v1_random	KI270717.1
chr14_KI270724v1_random	KI270724.1
chr9_KI270720v1_random	KI270720.1
chr14_KI270723v1_random	KI270723.1
chr9_KI270718v1_random	KI270718.1
chrUn_KI270317v1	KI270317.1
chrY_KI270740v1_random	KI270740.1
chrUn_KI270755v1	KI270755.1
chr1_KI270707v1_random	KI270707.1
chrUn_KI270579v1	KI270579.1
chrUn_KI270752v1	KI270752.1
chrUn_KI270512v1	KI270512.1
chrUn_KI270322v1	KI270322.1
chrUn_GL000226v1	GL000226.1
chrUn_KI270311v1	KI270311.1
chrUn_KI270366v1	KI270366.1
chrUn_KI270511v1	KI270511.1
chrUn_KI270448v1	KI270448.1
chrUn_KI270521v1	KI270521.1
chrUn_KI270581v1	KI270581.1
chrUn_KI270582v1	KI270582.1
chrUn_KI270515v1	KI270515.1
chrUn_KI270588v1	KI270588.1
chrUn_KI270591v1	KI270591.1
chrUn_KI270522v1	KI270522.1
chrUn_KI270507v1	KI270507.1
chrUn_KI270590v1	KI270590.1
chrUn_KI270584v1	KI270584.1
chrUn_KI270320v1	KI270320.1
chrUn_KI270382v1	KI270382.1
chrUn_KI270468v1	KI270468.1
chrUn_KI270467v1	KI270467.1
chrUn_KI270362v1	KI270362.1
chrUn_KI270517v1	KI270517.1
chrUn_KI270593v1	KI270593.1
chrUn_KI270528v1	KI270528.1
chrUn_KI270587v1	KI270587.1
chrUn_KI270364v1	KI270364.1
chrUn_KI270371v1	KI270371.1
chrUn_KI270333v1	KI270333.1
chrUn_KI270374v1	KI270374.1
chrUn_KI270411v1	KI270411.1
chrUn_KI270414v1	KI270414.1
chrUn_KI270510v1	KI270510.1
chrUn_KI270390v1	KI270390.1
chrUn_KI270375v1	KI270375.1
chrUn_KI270420v1	KI270420.1
chrUn_KI270509v1	KI270509.1
chrUn_KI270315v1	KI270315.1
chrUn_KI270302v1	KI270302.1
chrUn_KI270518v1	KI270518.1
chrUn_KI270530v1	KI270530.1
chrUn_KI270304v1	KI270304.1
chrUn_KI270418v1	KI270418.1
chrUn_KI270424v1	KI270424.1
chrUn_KI270417v1	KI270417.1
chrUn_KI270508v1	KI270508.1
chrUn_KI270303v1	KI270303.1
chrUn_KI270381v1	KI270381.1
chrUn_KI270529v1	KI270529.1
chrUn_KI270425v1	KI270425.1
chrUn_KI270396v1	KI270396.1
chrUn_KI270363v1	KI270363.1
chrUn_KI270386v1	KI270386.1
chrUn_KI270465v1	KI270465.1
chrUn_KI270383v1	KI270383.1
chrUn_KI270384v1	KI270384.1
chrUn_KI270330v1	KI270330.1
chrUn_KI270372v1	KI270372.1
chrUn_KI270548v1	KI270548.1
chrUn_KI270580v1	KI270580.1
chrUn_KI270387v1	KI270387.1
chrUn_KI270391v1	KI270391.1
chrUn_KI270305v1	KI270305.1
chrUn_KI270373v1	KI270373.1
chrUn_KI270422v1	KI270422.1
chrUn_KI270316v1	KI270316.1
chrUn_KI270338v1	KI270338.1
chrUn_KI270340v1	KI270340.1
chrUn_KI270583v1	KI270583.1
chrUn_KI270334v1	KI270334.1
chrUn_KI270429v1	KI270429.1
chrUn_KI270393v1	KI270393.1
chrUn_KI270516v1	KI270516.1
chrUn_KI270389v1	KI270389.1
chrUn_KI270466v1	KI270466.1
chrUn_KI270388v1	KI270388.1
chrUn_KI270544v1	KI270544.1
chrUn_KI270310v1	KI270310.1
chrUn_KI270412v1	KI270412.1
chrUn_KI270395v1	KI270395.1
chrUn_KI270376v1	KI270376.1
chrUn_KI270337v1	KI270337.1
chrUn_KI270335v1	KI270335.1
chrUn_KI270378v1	KI270378.1
chrUn_KI270379v1	KI270379.1
chrUn_KI270329v1	KI270329.1
chrUn_KI270419v1	KI270419.1
chrUn_KI270336v1	KI270336.1
chrUn_KI270312v1	KI270312.1
chrUn_KI270539v1	KI270539.1
chrUn_KI270385v1	KI270385.1
chrUn_KI270423v1	KI270423.1
chrUn_KI270392v1	KI270392.1
chrUn_KI270394v1	KI270394.1
__END__
