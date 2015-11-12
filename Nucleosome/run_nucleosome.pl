#!/usr/bin/perl

# usage: run_nucleosome.pl <fasta> <conc> <temp>

use strict;
use File::Basename;
my $start = time;


# chromosome fasta to work on is provided as first argument, e.g. chr1.fasta
my $file = shift @ARGV;
my $conc = shift @ARGV || undef; # their default value is 0.1
my $temp = shift @ARGV || undef; # their default value is 1

print "# Processing $file with concentration $conc and temperature $temp #\n";

# open files
my ($fname, $fpath, $fsuffix) = fileparse($file, qw(\.fa \.fasta \.fa\.gz \.fasta\.gz));

print "  Preparing input file $file...\n";
if ($fsuffix =~ /gz$/) {
	open(IN, "gzip -dc $file |") or die "   unable to open input $file\n";
}
else {
	open(IN, "$file") or die "   unable to open input $file\n";
}
my $fixfile = "$fname.fa";
print "  Opening fixed file $fixfile....\n";
open(OUT, ">$fixfile") or 
	die "unable to open output $fixfile\n";


# convert all Ns to As, skipping the sequence definition line
while (<IN>) {
	unless (/^>/) {
		s/N/A/gi;
	}
	print OUT "$_";
}

close IN;
close OUT;
printf "  Stripped $file in %.01f minutes\n", (time - $start)/60;



# run the nucleosome prediction program here
# generate title/name from the basename of the fasta file
my $name = $fname;
$name .= "_c$conc" if defined $conc;
$name .= "_t$temp" if defined $temp;

my @command = (
	# hard code paths, nucleosome in current directory
	# '/uufs/chpc.utah.edu/common/home/hcibcore/u0462865/perl/perls/5.16/bin/perl',
	'./nucleosome_prediction_v3/nucleosome_prediction.pl',
	'-t', "$name\_prediction",
	'-s', "$fixfile", 
	'-p', "$name\_prediction",
	'-tab',
);
if ($conc) {
	push @command, '-c', $conc;
}
if ($temp) {
	push @command, '-temp', $temp;
}

printf "  Executing command: %s\n", join(' ', @command);
system(@command) == 0 or die "   FAILED to execute command!\n";
	



# clean up the fasta files I made
unlink($fixfile, $file); 




# clean up the output file
# unfortunately, the script converts the fasta file to all caps before scanning,
# but this also includes the sequence names - doh!
# so chr1 becomes CHR1
print "  Cleaning up the output file....\n";
unless (-e "$name\_prediction.tab") {
	warn "Cannot find output file $name\_prediction.tab. Exiting.\n";
	exit;
}
my $infile = "$name\_prediction.tab";
my $outfile = "$name\_prediction.txt";

open(IN, $infile) or die "unable to open input file $infile!\n";
open(OUT, ">$outfile") or die "unable to open output file $outfile\n";

while (<IN>) {
	s/CHR/chr/;
	s/ZV9/Zv9/;
	s/SCAFFOLD/scaffold/;
	s/RANDOM/random/;
	s/UN/Un/;
	s/ULTRACONTIG/ultracontig/;
	print OUT "$_";
}
close IN;
close OUT;

# clean up input file
unlink $infile;

printf "  Finished with $name in %.01f minutes\n", (time - $start)/60;
