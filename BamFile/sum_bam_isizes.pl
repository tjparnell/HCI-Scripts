#!/usr/bin/perl

# a script to count paired-end insert sizes

use strict;
use List::Util qw(min max sum0);
use Bio::ToolBox::Data;
use Bio::ToolBox::db_helper 1.60 qw(
	open_db_connection 
	low_level_bam_fetch
);
use Bio::ToolBox::utility;


unless (@ARGV) {
	print " $0 <outfile> <bamfile1> \n";
	exit;
}

my $outfile = shift @ARGV;
my $file = shift @ARGV;

# count hashes
my %Fsize2count;
my %Rsize2count;

# processing callback
my $callback = sub {
	my $a = shift;
	
	# check alignment
	return if $a->unmapped;
	return unless $a->proper_pair;
	
	# record the insert size
	my $isize = $a->isize;
	if ($a->reversed) {
		$Rsize2count{$isize} += 1;
	}
	else {
		$Fsize2count{$isize} += 1;
	}
};


print "counting $file...\n";
my $sam = open_db_connection($file) or
	die "can't open bam file $file\n";

# Loop through the chromosomes
for my $tid (0 .. $sam->n_targets - 1) {
	low_level_bam_fetch($sam, $tid, 0, $sam->target_len($tid), $callback, 1);
}
	
# now record the size counts 
my $minsize = min(keys(%Fsize2count), keys(%Rsize2count));
my $maxsize = max(keys(%Fsize2count), keys(%Rsize2count));
my $Data = Bio::ToolBox::Data->new(
	feature => 'alignment_summary', 
	columns => ['size', 'F_isize', 'R_isize'],
) or die "unable to generate Data structure!";
for my $size ($minsize .. $maxsize) {
	my $f = $Fsize2count{$size} || 0;
	my $r = $Rsize2count{$size} || 0;
	next if ($f == 0 and $r == 0); # skip empties
	$Data->add_row([$size, $f, $r]);
}

# Finish
my $success_write = $Data->save(
	'filename' => $outfile,
);
if ($success_write) {
	print "wrote $success_write!\n";
}

# print summary
printf "%12s Forward alignments had a positive insert size\n",
    format_with_commas(sum0( map {$Fsize2count{$_}} grep {$_ > 0} keys %Fsize2count));
printf "%12s Forward alignments had a negative insert size\n",
    format_with_commas(sum0( map {$Fsize2count{$_}} grep {$_ < 0} keys %Fsize2count));
printf "%12s Reverse alignments had a positive insert size\n",
    format_with_commas(sum0( map {$Rsize2count{$_}} grep {$_ > 0} keys %Rsize2count));
printf "%12s Reverse alignments had a negative insert size\n",
    format_with_commas(sum0( map {$Rsize2count{$_}} grep {$_ < 0} keys %Rsize2count));




