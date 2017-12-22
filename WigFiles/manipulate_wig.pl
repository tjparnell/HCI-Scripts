#!/usr/bin/env perl

# Timothy J. Parnell, PhD
# Huntsman Cancer Institute
# University of Utah
# Salt Lake City, UT 84112
#  
# This package is free software; you can redistribute it and/or modify
# it under the terms of the Artistic License 2.0.  
# 
# Updated versions of this file may be found in the repository
# https://github.com/tjparnell/HCI-Scripts/

use strict;
use Getopt::Long;
use IO::File;
use IO::Handle;

# version
# 1.0 initial release
# 1.1 add absolute option and optimize parsing
# 1.2 streamline wig format processing, option to skip chromosomes
# 1.3 bug fixes
# 1.4 add chromosome apply regex
my $VERSION = 1.4;



unless (@ARGV) {
	print <<END;

A script to manipulate the score value of wig files. This will process all 
forms of text based wig files, including fixedStep, variableStep, and bedGraph. 
Files can be gzipped. 

NOTE: More than one option may be specified! The options below are the order 
in which the score is manipulated. If they are not in the order you want, you 
may have to pipe to sequential instances. Use 'stdin' and 'stdout' for filenames.
Use an equal sign to define options with negative values, e.g. --mult=-1

BIGWIG files: The UCSC bigWigToWig and wigToBigWig utilities accept 'stdin' 
and 'stdout' as file names, so you can pipe out and in to bigWig formats. 

Usage: $0 [options] -i <file1.wig> -o <file1.out.wig>
Options: 
  --in <file>      Input file. Accepts gz compression. Accepts 'stdin'.
  --out <file>     Output file. Accepts gz compression. Accepts 'stdout'.
                   Optional if all you want is to calculate statistics.
  --skip <regex>   Discard lines where chromosomes match the regular 
                   expression. Example: 'chrM|chrUn|random'
                   Good for removing chromosomes from analysis and/or wig.
  --apply <regex>  Apply manipulations or statistics ONLY to chromosomes 
                   that match the regular expression, leaving other lines 
                   untouched. Example: 'chrX'
  --null           Convert null, NA, N/A, NaN, inf values to 0
  --delog <int>    Delog values of base [int], usually 2 or 10
  --abs            Convert to the absolute value 
  --mult <float>   Multiply score by the given value
  --add <float>    Add the given value to the score
  --log [2|10]     Convert to log2 or log10. Scores of 0 are left as 0.
  --min <float>    Set the minimum score
  --max <float>    Set the maximum score
  --place <int>    Format the score to the given number of decimal positions
  --zero           Discard lines with zero values
  --stats          Calculate statistics across the genome at base pair resolution.
                   Statistics are calculated after the above processing. Only 
                   applied chromosomes are calculated.
END
	exit;
}

# Options
my $infile;
my $outfile;
my $skip;
my $apply;
my $doNull = 0;
my $deLogValue;
my $doAbsolute = 0;
my $multiplyValue;
my $addValue;
my $logValue;
my $places;
my $minValue;
my $maxValue;
my $noZeroes;
my $doStats;
GetOptions( 
	'input=s'       => \$infile,
	'output=s'      => \$outfile,
	'skip=s'        => \$skip,
	'apply=s'       => \$apply,
	'null!'         => \$doNull,
	'delog=i'       => \$deLogValue,
	'abs!'          => \$doAbsolute,
	'multiply=f'    => \$multiplyValue,
	'add=f'         => \$addValue,
	'log=i'         => \$logValue,
	'place=i'       => \$places,
	'minimum=f'     => \$minValue, # 
	'maximum=f'     => \$maxValue, # 
	'zero'          => \$noZeroes,
	'stats!'        => \$doStats,
) or die "bad options!\n";


# checks
die "no input file provided!\n" unless $infile;
my $doMin = defined $minValue ? 1 : 0;
my $doMax = defined $maxValue ? 1 : 0;
if ($logValue) {
	$logValue = $logValue == 2 ? log(2) : $logValue == 10 ? log(10) : undef;
	die "bad log value!\n" unless defined $logValue;
}
if (defined $places) {
	$places = '%.' . $places . 'f';
}

# chromosome skipping regex
my ($skip_regex, $apply_regex);
if ($skip) {
	$skip_regex = qr($skip);
}
if ($apply) {
	$apply_regex = qr($apply);
}


# open file handles
my ($infh, $outfh);
if ($infile =~ /^stdin$/i) {
	$infh = IO::Handle->new;
	$infh->fdopen(fileno(STDIN), 'r');
}
elsif (-e $infile) {
	$infh = $infile =~ /\.gz$/i ? 
			IO::File->new("gzip -dc $infile |") :
			IO::File->new($infile);
	die "can't open $infile! $!" unless $infh;
}
else {
	die "unrecognized $infile!";
}
if ($outfile =~ /^stdout$/i) {
	$outfh = IO::Handle->new;
	$outfh->fdopen(fileno(STDOUT), 'w');
}
elsif ($outfile) {
	$outfh = $outfile =~ /\.gz$/i ? 
			IO::File->new("| gzip >>$outfile") :
			IO::File->new($outfile, 'w');
	die "can't open $outfile! $!" unless $outfh;
}


# stats hash
my $stats = {
	count         => 0,
	sumData       => 0,
	sumSquares    => 0,
	minVal        => undef,
	maxVal        => undef,
};


# walk through the file
my $count = 0;
my $span = 1;
my $chrom_skip = 0;
my $chrom_ignore = 0;
my $wig_process_sub;
while (my $line = $infh->getline) {
	# look at the first characters to determine the type of line we have
	my $prefix = lc substr($line,0,5);
	if ($prefix eq 'track') {
		# track line
		$outfh->print($line) if $outfh;
		next;
	}
	elsif ($prefix eq 'brows') {
		# browser line
		$outfh->print($line) if $outfh;
		next;
	}
	elsif ($prefix eq 'varia' or $prefix eq 'fixed') {
		# a step definition line
		if ($line =~ /chrom=([\w\-\.]+)/) {
			# check the chromosome
			my $chrom = $1;
			if ($skip_regex and $chrom =~ $skip_regex) {
				print STDERR "skipping chromosome $chrom\n";
				$chrom_skip = 1;
			}
			else {
				$chrom_skip = 0;
			}
			if ($apply and $chrom !~ $apply_regex) {
				print STDERR "ignoring chromosome $chrom\n";
				$chrom_ignore = 1;
			}
			else {
				$chrom_ignore = 0;
			}
		}
		if ($line =~ /span=(\d+)/i) {
			# capture span size if present
			$span = $1;
		}
		$outfh->print($line) if $outfh;
		next;
	} 
	elsif (substr($prefix,0,1) eq '#') {
		# comment line
		$outfh->print($line) if $outfh;
		next;
	}
	
	# skipping current chromosome
	next if $chrom_skip;
	
	# ignoring current chromosome
	if ($chrom_ignore) {
		$outfh->print($line) if $outfh;
		next;
	}
	
	# determine format
	unless (defined $wig_process_sub) {
		my @data = split /\s+/, $line;
		if (scalar @data == 4) {
			$wig_process_sub = \&process_bedGraph;
		}
		elsif (scalar @data == 2) {
			$wig_process_sub = \&process_variableStep;
		}
		elsif (scalar @data == 1) {
			$wig_process_sub = \&process_fixedStep;
		}
	}
	
	# process
	chomp $line;
	&$wig_process_sub($line);
}


# close up shop
$infh->close;
$outfh->close if $outfh;

# print final messages
my $statMessage;
if ($doStats) {
	my $basecount = $stats->{count};
	my $min   = $stats->{minVal};
	my $max   = $stats->{maxVal};
	my $mean  = $stats->{count} ? sprintf("%.05f", $stats->{sumData} / $stats->{count}) : 0;
	my $stddev = sprintf("%.05f", sqrt(binVariance()) );
	$statMessage = <<STATS; 
basesCovered: $basecount
mean: $mean
min: $min
max: $max
std: $stddev
STATS
}

if ($outfile =~ /stdout/i) {
	print STDERR " converted $count lines, wrote file $outfile\n" if $outfile;
	print STDERR $statMessage if $statMessage;
}
else {
	print STDOUT " converted $count lines, wrote file $outfile\n" if $outfile;
	print STDOUT $statMessage if $statMessage;
}

sub process_bedGraph {
	my @data = split "\t", shift;
	return if ($skip and $data[0] =~ $skip_regex);
	if ($apply and $data[0] !~ $apply_regex) {
		$outfh->printf("%s\t%d\t%d\t%s\n", $data[0], $data[1], $data[2], $data[3]) 
			if $outfh;
		return;
	}
	$data[3] = process_score($data[3]);
	return if not defined $data[3];
	$outfh->printf("%s\t%d\t%d\t%s\n", $data[0], $data[1], $data[2], $data[3]) 
		if (defined $data[3] and $outfh);
	$count++;
	if ($doStats and defined $data[3]) {
		my $length = $data[2] - $data[1];
		$stats->{count} += $length;
		$stats->{sumData} += ($length * $data[3]);
		$stats->{sumSquares} += ( ($data[3] ** 2) * $length );
		$stats->{minVal} = $data[3] if not defined $stats->{minVal};
		$stats->{minVal} = $data[3] if $data[3] < $stats->{minVal};
		$stats->{maxVal} = $data[3] if not defined $stats->{maxVal};
		$stats->{maxVal} = $data[3] if $data[3] > $stats->{maxVal};
	}
}

sub process_variableStep {
	my @data = split /\s+/, shift; # could be either tab or space
	$data[1] = process_score($data[1]);
	$outfh->printf("%d %s\n", $data[0], $data[1]) if (defined $data[1] and $outfh);
	$count++;
	process_step_stats($data[1]) if $doStats;
}

sub process_fixedStep {
	my $score = shift;
	$score = process_score($score);
	$outfh->printf("%s\n",$score) if (defined $score and $outfh);
	$count++;
	process_step_stats($score) if $doStats;
}

sub process_score {
	my $v = shift; # score
	if ($doNull and $v =~ /^(?:n.?[na])|(?:\-?inf)/i) {$v = 0}
	if ($deLogValue) {$v = $deLogValue ** $v}
	if ($doAbsolute) {$v = abs($v)}
	if ($multiplyValue) {$v *= $multiplyValue}
	if ($addValue) {$v += $addValue}
	if ($logValue) {$v = $v == 0 ? 0 : log($v) / $logValue}
	if ($doMin and $v < $minValue) {$v = $minValue}
	if ($doMax and $v > $maxValue) {$v = $maxValue}
	if ($places) {$v = sprintf($places, $v)};
	return undef if ($noZeroes and $v == 0);
	return $v;
}

sub process_step_stats {
	return unless defined $_[0];
	for (1 .. $span) {
		$stats->{count} += 1;
		$stats->{sumData} += $_[0];
		$stats->{sumSquares} += $_[0] ** 2;
		$stats->{minVal} = $_[0] if not defined $stats->{minVal};
		$stats->{maxVal} = $_[0] if not defined $stats->{maxVal};
		$stats->{minVal} = $_[0] if $_[0] < $stats->{minVal};
		$stats->{maxVal} = $_[0] if $_[0] > $stats->{maxVal};
	}
}

sub binVariance {
    return 0 unless $stats->{count};
    my $var = $stats->{sumSquares} - $stats->{sumData}**2/$stats->{count};
    if ($stats->{count} > 1) {
	$var /= $stats->{count}-1;
    }
    return 0 if $var < 0;
    return $var;
}



