#!/usr/bin/perl

use strict;
use Bio::ToolBox::Data;



unless (@ARGV) {
	print <<DOC;
A quick script to compare two VCF files with each other.
Looks for identical variants and reports the intersection.
This is done for named variants (e.g. dbSNP names), and novel 
variants (no name). It will also report the numbers of named, 
novel, and filtered variants in each file.

$0 <file1.vcf> <file2.vcf> <outfile.txt>

DOC
	exit;
}


my ($vcf1, $vcf2, $outfile) = @ARGV or die "need three files!\n";

# variant hashes
# {seq_id}{start}{alt} = 1, 2, or 3 (both)
my %named;
my %unnamed;
my %filtered;



### first file
my $Stream1 = Bio::ToolBox::Data->new(
	file   => $vcf1,
	stream => 1,
) or die "unable to open $vcf1! $!\n";
unless ($Stream1->vcf) {
	die "$vcf1 doesn't appear to a be a vcf file!\n";
}

while (my $row = $Stream1->next_row) {
	my $s = $row->seq_id;
	my $p = $row->start;
	my $a = $row->value(4); # alternate
	my $f = $row->value(6); # filter
	if ($f eq '.' or $f eq 'PASS') {
		# we have an assumed or explicit passed variant
		if ($row->value(2) eq '.') {
			# an unnamed variant
			$unnamed{$s}{$p}{$a} = 1;
		}
		else {
			# a named variant
			# we're not going to check names, since we're assuming coordinate and alt 
			# is good enough identifiers here, and the names in theory should be identical
			$named{$s}{$p}{$a} = 1;
		}
	} else {
		# we have a failed variant
		$filtered{$s}{$p}{$a} = 1;
	}
}
$Stream1->close_fh;



### second file
my $Stream2 = Bio::ToolBox::Data->new(
	file   => $vcf2,
	stream => 1,
) or die "unable to open $vcf2! $!\n";
unless ($Stream2->vcf) {
	die "$vcf2 doesn't appear to a be a vcf file!\n";
}

while (my $row = $Stream2->next_row) {
	my $s = $row->seq_id;
	my $p = $row->start;
	my $a = $row->value(4); # alternate
	my $f = $row->value(6); # filter
	if ($f eq '.' or $f eq 'PASS') {
		# we have an assumed or explicit passed variant
		if ($row->value(2) eq '.') {
			# an unnamed variant
			# check for existing
			if (exists $unnamed{$s}{$p}{$a}) {
				$unnamed{$s}{$p}{$a} = 3;
			}
			else {
				$unnamed{$s}{$p}{$a} = 2;
			}
		}
		else {
			# a named variant
			# we're not going to check names, since we're assuming coordinate and alt 
			# is good enough identifiers here, and the names in theory should be identical
			# check for existing
			if (exists $named{$s}{$p}{$a}) {
				$named{$s}{$p}{$a} = 3;
			}
			else {
				$named{$s}{$p}{$a} = 2;
			}
		}
	} else {
		# we have a failed variant
		# check for existing
		if (exists $filtered{$s}{$p}{$a}) {
			$filtered{$s}{$p}{$a} = 3;
		}
		else {
			$filtered{$s}{$p}{$a} = 2;
		}
	}
}
$Stream2->close_fh;



### Compile results
my $Data = Bio::ToolBox::Data->new(
	columns => [qw(File Total Named Novel Filtered Named_Unique Named_Shared 
					Novel_Unique Novel_Shared Filtered_Unique Filtered_Shared)]
);
my @stats1 = ($vcf1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
my @stats2 = ($vcf2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);

# count up named variants
foreach my $s (keys %named) {
	foreach my $p (keys %{$named{$s}}) {
		foreach my $a (keys %{ $named{$s}{$p} }) {
			my $v = $named{$s}{$p}{$a};
			if ($v == 1) {
				$stats1[1]++; # Total
				$stats1[2]++; # Named
				$stats1[5]++; # Named_Unique
			} elsif ($v == 2) {
				$stats2[1]++; # Total
				$stats2[2]++; # Named
				$stats2[5]++; # Named_Unique
			} else {
				$stats1[1]++; # Total
				$stats2[1]++; # Total
				$stats1[2]++; # Named
				$stats2[2]++; # Named
				$stats1[6]++; # Named_Shared
				$stats2[6]++; # Named_Shared
			}
		}
	}
}

# count up novel variants
foreach my $s (keys %unnamed) {
	foreach my $p (keys %{$unnamed{$s}}) {
		foreach my $a (keys %{ $unnamed{$s}{$p} }) {
			my $v = $unnamed{$s}{$p}{$a};
			if ($v == 1) {
				$stats1[1]++; # Total
				$stats1[3]++; # Novel
				$stats1[7]++; # Novel_Unique
			} elsif ($v == 2) {
				$stats2[1]++; # Total
				$stats2[3]++; # Novel
				$stats2[7]++; # Novel_Unique
			} else {
				$stats1[1]++; # Total
				$stats2[1]++; # Total
				$stats1[3]++; # Novel
				$stats2[3]++; # Novel
				$stats1[8]++; # Novel_Shared
				$stats2[8]++; # Novel_Shared
			}
		}
	}
}

# count up filtered variants
foreach my $s (keys %filtered) {
	foreach my $p (keys %{$filtered{$s}}) {
		foreach my $a (keys %{ $filtered{$s}{$p} }) {
			my $v = $filtered{$s}{$p}{$a};
			if ($v == 1) {
				$stats1[1]++; # Total
				$stats1[4]++; # Filtered
				$stats1[9]++; # Filtered_Unique
			} elsif ($v == 2) {
				$stats2[1]++; # Total
				$stats2[4]++; # Filtered
				$stats2[9]++; # Filtered_Unique
			} else {
				$stats1[1]++; # Total
				$stats2[1]++; # Total
				$stats1[4]++; # Filtered
				$stats2[4]++; # Filtered
				$stats1[10]++; # Filtered_Shared
				$stats2[10]++; # Filtered_Shared
			}
		}
	}
}

# store
$Data->add_row(\@stats1);
$Data->add_row(\@stats2);
my $success = $Data->save($outfile);
print "wrote $success\n";


