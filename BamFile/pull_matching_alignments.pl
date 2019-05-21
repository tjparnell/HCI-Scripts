#!/usr/bin/env perl

use strict;
use Bio::ToolBox 1.65;
use Bio::ToolBox::db_helper;

unless (@ARGV) {
	print <<USAGE;

Extract sam alignments based on name.

Usage: $0 <input_bam> <list> <output_bam>

  <input_bam> is the source bam file
  <list> is a text file containing a list of read (query) names, one name per line.
  <output_bam> is the output bam file. This is optional

Alignments from the input bam file with a matching name are written to a new bam file.
USAGE
	exit 0;
}

# arguments
my $infile   = shift @ARGV;
my $namefile = shift @ARGV;
my $outfile  = shift @ARGV;

# check for required
die "need a list file!\n" unless $namefile;
unless ($outfile) {
	$outfile = $infile;
	$outfile =~ s/\.bam$/_matched.bam/i;
}


# load up list file
my %lookup;
my $duplicates = 0;
my $fh = Bio::ToolBox->open_file($namefile) or 
	die "unable to open file! $!\n";
while (my $line = $fh->getline) {
	chomp $line;
	if (exists $lookup{$line}) {
		$duplicates++;
	}
	else {
		$lookup{$line} = 0;
	}
}
$fh->close;
print " Warning: there were $duplicates duplicate query names in the list file!\n" if $duplicates;
printf "Loaded %s query names\n", scalar keys %lookup;

# open input bam file
my $sam = Bio::ToolBox->open_database($infile) or 
	die " Cannot open input Bam file!\n";

# read header and set adapter specific alignment writer
# also get low level bam adapter
my ($lowbam, $header, $write_alignment);
if (Bio::ToolBox->bam_adapter eq 'sam') {
	$lowbam = $sam->bam;
	$header = $lowbam->header;
	$write_alignment = \&write_sam_alignment;
}
elsif (Bio::ToolBox->bam_adapter eq 'hts') {
	$lowbam = $sam->hts_file;
	$header = $lowbam->header_read;
	$write_alignment = \&write_hts_alignment;
}
else {
	die "unrecognized bam adapter!"; # shouldn't happen
}

# open output bam file
my $out = Bio::ToolBox::db_helper::write_new_bam_file($outfile) or 
	die "cannot open output bam file '$outfile'!\n";
$out->header_write($header);




# look for alignments
my $matched = 0;
while (my $a = $lowbam->read1($header)) {
	if (exists $lookup{ $a->qname }) {
		$matched++;
		$lookup{ $a->qname }++;
		&$write_alignment($a);
	}
}

# check for leftovers
my $notfound = 0;
foreach (keys %lookup) {
	$notfound++ if $lookup{$_} == 0;
}
print "Warning: $notfound names could not be found!\n" if $notfound;

print "Matched $matched alignments\n";
print "Wrote output file $outfile\n";
exit 0;


sub write_sam_alignment {
	# wrapper for writing Bio::DB::Sam alignments
	return $out->write1($_[0]);
}

sub write_hts_alignment {
	# wrapper for writing Bio::DB::HTS alignments
	return $out->write1($header, $_[0]);
}


