#!/usr/bin/perl

# a script to count paired-end insert sizes

use strict;

use Bio::ToolBox::Data;
use Bio::ToolBox::db_helper::bam;


# variables
# these could be obtained from a commandline if I had the wherewithall

my $min = 1;
my $max = 500;

unless (@ARGV) {
	die " $0 <outfile> <bamfile1> <bamfile2>....\n";
}

my $outfile = shift @ARGV;
my @files = @ARGV;

# make output structure
my $Data = Bio::ToolBox::Data->new(
	feature => 'alignment_summary', 
	columns => ['size', @files],
) or die "unable to generate Data structure!";
my %size2count;

# fill up Data structures with sizes
for (my $size = $min; $size <= $max; $size++) {
	my $r = $Data->add_row;
	$Data->value($r, 0, $size);
	$size2count{$size} = 0;
}

# processing callback
my $callback = sub {
	my $a = shift;
	
	# check alignment
	return if $a->unmapped;
	return unless $a->proper_pair;
	
	# we only need to process one of the two pairs, 
	# so only take the left (forward strand) read
	return unless $a->strand == 1;
	
	# record the insert size
	my $isize = $a->isize;
	$size2count{$isize} += 1;
};

# collect Data from files
my $index = 1;
foreach my $file (@files) {
	print "counting $file...\n";
	my $sam = open_bam_db($file) or
		die "can't open bam file $file\n";
	
	# Loop through the chromosomes
	for my $tid (0 .. $sam->n_targets - 1) {
		
		# sequence name
		my $seq_id = $sam->target_name($tid);
		
		# process the reads
		$sam->fetch($seq_id, $callback);
	}
	
	# now record the size counts 
	$Data->iterate( sub {
		my $row = shift;
		my $size = $row->value(0);
		$row->value($index, $size2count{$size});
		# erase the count for next
		$size2count{$size} = 0;
	} );
	
	# prepare for next file
	$index++;
}

my $success_write = $Data->save(
	'filename' => $outfile,
);
if ($success_write) {
	print "wrote $success_write!";
}





