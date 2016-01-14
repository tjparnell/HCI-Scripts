#!/usr/bin/perl

use strict;
use Getopt::Long;
use Statistics::Lite qw(sum);
use Bio::ToolBox::Data;
use Bio::ToolBox::utility;


print "Sums the counts for individual and repeat_Family members\n";

unless (@ARGV) {
print <<USAGE ;
$0 -n <i> -f <i> -d <i> -o outfile -i infile

	-n i    0-based index for the Name column
	-f i    0-based index for the FamilyName column, as "Type_Family"
	-d i    0-based index for one or more columns of read counts to sum up
	        can give as comma-delimited list or range, ex 5,7-21
	-i file The input file name
	-o file The output file name 
	
USAGE
exit;
}

my ($name_i, $fam_i, $dlist, $infile, $outfile);
GetOptions(
	'n=i'   => \$name_i,
	'f=i'   => \$fam_i,
	'd=s'   => \$dlist,
	'i=s'   => \$infile,
	'o=s'   => \$outfile,
) or die "something with options is not right!";

my @list = parse_list($dlist);
unless (defined $name_i and defined $fam_i and scalar @list > 0) {
	die "Need at least three or more columns, fool!";
}

my $Stream = Bio::ToolBox::Data->new(
	in      => $infile,
	stream  => 1,
) or die "unable to open Data stream! $!\n";

# collect the data
my %data;
while (my $row = $Stream->next_row) {
	my $name = $row->value($name_i);
	$name =~ s/\.\d+$//;
	my ($type, $family) = split(/_/, $row->value($fam_i));
	my $s = sum( map {$row->value($_)} @list);
	
	# record
	$data{$type}{$family}{$name} += $s;
	$data{$type}{$family}{' ALL'} += $s;
	$data{$type}{' ALL'} += $s;
	$data{' ALL'} += $s;
}

# output the data
my $D = Bio::ToolBox::Data->new(
	columns => ['Type', 'Family', 'Name', 'Sum_' . join('_', map {$Stream->name($_)} @list) ],
);


foreach my $type (sort {$a cmp $b} keys %data) {
	if ($type eq ' ALL') {
		$D->add_row( ['All', 'ALL', 'ALL', $data{$type}] );
		next;
	}
	foreach my $family (sort {$a cmp $b} keys %{$data{$type}}) {
		if ($family eq ' ALL') {
			$D->add_row( [$type, 'ALL', 'ALL', $data{$type}{$family}] );
			next;
		}
		foreach my $name (sort {$a cmp $b} keys %{$data{$type}{$family}}) {
			if ($name eq ' ALL') {
				$D->add_row( [$type, $family, 'ALL', $data{$type}{$family}{$name}] );
				next;
			}
			$D->add_row( [ $type, $family, $name, $data{$type}{$family}{$name}] );
		}
	}
}

my $success = $D->save($outfile);
print "wrote file $success\n" if $success;
$Stream->close_fh;


