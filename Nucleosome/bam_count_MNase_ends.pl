#!/usr/bin/perl

use strict;
use Bio::ToolBox::Data;
use Bio::ToolBox::db_helper::bam;
use Bio::DB::Sam; # for fai index
my $VERSION = 'beta';

my $description = <<END;
calculate the cut site frequency for MNaseSeq
  grabs the 2 bp genomic sequence representing the -1 and 1 bp of each alignment
  reports the occurency of each
  MNase should cut predominantly at AA/AT/TA/TT
Usage: $0 <bamfile.bam> <genome.fa> <outfile.txt>
END

unless (@ARGV) {
	print $description;
	exit;
}


# parameters
my $bamfile = shift @ARGV or die "no bam file provided!\n $description\n";
my $fafile = shift @ARGV or die "no fasta file provided!\n $description\n";
my $outfile = shift @ARGV or die "no output file provided!\n $description\n";

# open files
my $sam = open_bam_db($bamfile) or die 
	" unable to open bam file '$bamfile'!\n";
my $fai = Bio::DB::Sam::Fai->open($fafile) or die 
	" unable to open fasta file '$fafile'!\n";

# loop through alignments on each chromosome
my %seq2count; 
for my $tid (0 .. $sam->n_targets - 1) {
	my $seq_id = $sam->target_name($tid);
	my $seq_length = $sam->target_len($tid);
	$seq2count{seqid} = $seq_id;
	
	print "  checking $seq_id....\n";
	$sam->bam_index->fetch($sam->bam, $tid, 0, $seq_length, \&callback, \%seq2count);
}
delete $seq2count{seqid}; # we don't need this anymore

# count up total
my $total = 0;
foreach (values %seq2count) { $total += $_}

# create output data and store counts
my $Data = Bio::ToolBox::Data->new(columns => [qw(Sequence Count Fraction)]);
foreach my $seq (sort {$a cmp $b} keys %seq2count) {
	$Data->add_row( [$seq, $seq2count{$seq}, sprintf("%.4f", $seq2count{$seq}/$total) ] );
}
my $success = $Data->save($outfile);
print "wrote file $success\n";	

# alignment callback
sub callback {
	my ($a, $s2c) = @_;
	my $s = $a->reversed ? $a->calend : $a->pos; 
		# returned coordinates are interbase
		# doesn't matter for reverse alignments, but forward alignments is actually -1
		# that's fine, I don't have to add and subtract 1 :-)
	my $e = $s + 1;
	my $seq = $fai->seq($s2c->{seqid}, $s, $e);
	$s2c->{$seq} += 1;
}

