#!/usr/bin/perl

use strict;
use File::Spec;
use Bio::DB::Sam;
use Bio::ToolBox::Data;

unless (@ARGV) {
	print "A script to count alignments for each sequence in bam files\n";
	print "Bam files do NOT need to be sorted or indexed.\n";
	print "Usage: $0 <outfile> <file1.bam ...>\n";
	exit;
}
my $outfile = shift @ARGV;

# Prepare data structure
my $Data = Bio::ToolBox::Data->new(
	datasets => [ qw(Reference) ],
	feature  => 'Reference_sequence',
);


while (@ARGV) {
	my $path = shift @ARGV;
	my (undef, $dir, $file) = File::Spec->splitpath($path);
	print "working on $path...\n";
	
	unless ($file =~ /\.bam$/i) {
		warn "$path is not a bam file!\n";
		next;
	}
	
	my $bam = Bio::DB::Bam->open($path);
	my $header = $bam->header; # must always get before reading alignments
	my ($ref2count, $id2ref) = set_up_header_hash($header);
	
	# check reference list
	unless ($Data->last_row) {
		add_reference_list($Data, $ref2count);
	}
	
	# count alignments
	while (my $a = $bam->read1) {
		$ref2count->{'total'}++;
		if ($a->unmapped) {
			$ref2count->{'unmapped'}++;
			next;
		}
		my $ref = $id2ref->{ $a->tid };
		$ref2count->{$ref}++;
	}
	
	# record data
	$file =~ s/\.bam$//i;
	my $index = $Data->add_column($file);
	$Data->metadata($index, 'source', $path);
	$Data->metadata($index, 'total', $ref2count->{'total'});
	my $stream = $Data->row_stream;
	while (my $row = $stream->next_row) {
		my $ref = $row->value(0); # always index 0
		$row->value($index, $ref2count->{$ref});
	}
	$Data->save(filename => $outfile);
}

print "Finished\n";



sub set_up_header_hash {
	my $head = shift;
	
	my @header = split("\n", $head->text);
	my %ref2count;
	my %id2ref;
	my $tid = 0;
	foreach (@header) {
		next unless /^\@SQ/;
		my ($type, $name, $length) = split /\s+/;
		$name =~ s/^SN://;
		$ref2count{$name} = 0;
		$id2ref{$tid} = $name;
		$tid++;
	}
	$ref2count{'total'} = 0;
	$ref2count{'unmapped'} = 0;
	return (\%ref2count, \%id2ref);
}


sub add_reference_list {
	my ($Data, $ref2count) = @_;
	
	# add first lines
	$Data->add_row( [qw(unmapped)] );
	
	my @list = sort {$a cmp $b} keys %$ref2count;
	while (@list) {
		my $ref = shift @list;
		next if ($ref =~ /total|unmapped/i);
		$Data->add_row( [$ref] );
	}
}


