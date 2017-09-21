#!/usr/bin/perl

# dump feature

use strict;
use Getopt::Long;
use Data::Dumper;
use Bio::ToolBox::db_helper qw(
	get_db_feature
);

my $help = <<HELP;
Dump db features to stdout as text
Usage $0 
   --db database
   --name name
   --type type
   --id primary_id
HELP

unless (@ARGV) {
	print $help;
	exit;
}

my $database;
my $name;
my $type;
my $id;
GetOptions(
	'db=s'      => \$database,
	'name=s'    => \$name,
	'type=s'    => \$type,
	'id=i'      => \$id,
) or die $help;

my $f = get_db_feature(
	db   => $database,
	id   => $id,
	name => $name,
	type => $type,
);

if ($f) {
# 	my @s = $f->get_SeqFeatures;
# 	$f->add_segment(@s); # this will add the subfeatures, but is not recursive
	print Dumper($f) . "\n";
# 	foreach (@s) {print Dumper($_) . "\n"}
	# my $f2 = $f->clone;
	# print Dumper($f2) . "\n";
}
else {
	print " Feature not found\n";
}


