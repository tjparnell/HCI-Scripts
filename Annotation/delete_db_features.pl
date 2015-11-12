#!/usr/bin/env perl

use strict;
use Bio::DB::SeqFeature::Store;


print "\n a script to delete feature types from a Bio::DB::SeqFeature::Store database\n";
unless (@ARGV) {
	print "\n usage: $0 <database> <type> <type> ...\n";
	print "   WARNING: make sure you know what you are doing!!!\n";
	exit;
}

my $database = shift @ARGV;

my $db = Bio::DB::SeqFeature::Store->new(
	-adaptor => 'DBI::mysql',
	-dsn     => $database,
	-user    => 'loader',
	-pass    => '',
) or die "can't connect to $database!\n";


while (@ARGV) {
	my $type = shift @ARGV;
	print "\n deleting $type from $database...\n";
	
	my $stream = $db->get_seq_stream(
		-types   => $type,
	);
	
	my $deleted = 0;
	my $not_deleted = 0;
	while (my $f = $stream->next_seq) {
		my $success = $db->delete($f);
		if ($success) {
			$deleted++;
		}
		else {
			$not_deleted++;
		}
		if ($deleted % 1000 == 0) {
			print "   $deleted features deleted...\n";
		}
	}
	print " $deleted features were deleted\n" if $deleted;
	print " $not_deleted features were not deleted\n" if $not_deleted;
}