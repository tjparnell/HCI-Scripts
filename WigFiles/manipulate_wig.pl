#!/usr/bin/env perl

use strict;
use Getopt::Long;
use IO::File;
use IO::Handle;

unless (@ARGV) {
	print <<END;

A script to manipulate the score value of wig files. This will process all 
forms of text based wig files, including fixedStep, variableStep, and bedGraph. 
Files can be gzipped. 

NOTE: More than option may be specified! The options below are the order in 
which the score is manipulated. If they are not in the order you want, you 
may have to pipe to sequential instances. Use 'stdin' and 'stdout' for 
filenames.

Usage: $0 [options] -i <file1.wig> -o <file1.out.wig>
Options: 
  --in <file>      Input file. Accepts gz compression. Accepts 'stdin'.
  --out <file>     Output file. Accepts gz compression. Accepts 'stdout'.
                   Optional if all you want is to calculate statistics.
  --null           Convert null, NA, N/A, NaN, inf values to 0
  --delog <int>    Delog values of base [int], usually 2 or 10
  --flip           Flip the sign of the score value. All negative scores 
                   become positive, and vice versa. 
  --mult <float>   Multiply score by the given value
  --add <float>    Add the given value to the score
  --log [2|10]     Convert to log2 or log10. Scores of 0 are left as 0.
  --min <float>    Set the minimum score
  --max <float>    Set the maximum score
  --place <int>    Format the score to the given number of decimal positions
  --zero           Discard lines with zero values
  --stats          Calculate statistics across the genome at base pair resolution.
                   Statistics are calculated after above processing.
END
	exit;
}

# Options
my $infile;
my $outfile;
my $doNull = 0;
my $deLogValue;
my $doFlip = 0;
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
	'null!'         => \$doNull,
	'delog=i'       => \$deLogValue,
	'flip!'         => \$doFlip,
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
while (my $line = $infh->getline) {
	$line =~ s/[\r\n]+$//; # strip all line endings
	my @data = split /\s+/, $line;
	
	# a track line
	if ($data[0] =~ /^track|browser/i) {
		$outfh->print("$line\n") if $outfh;
	}
	
	# a step definition line
	elsif ($data[0] =~ /^variablestep|fixedstep/i) { 
		if ($line =~ /span=(\d+)/i) {
			# capture span size if present
			$span = $1;
		}
		$outfh->print("$line\n") if $outfh;
	} 
	
	# comment line
	elsif ($data[0] =~ /^#/) {
		$outfh->print("$line\n") if $outfh;
	}
	
	# a bedGraph data line
	elsif (scalar @data == 4) {
		$data[3] = process_score($data[3]);
		$outfh->print(join("\t", @data) . "\n") if (defined $data[3] and $outfh);
		$count++;
		process_bdg_stats(@data) if $doStats;
	} 
	
	# a variableStep data line
	elsif (scalar @data == 2) { 
		$data[1] = process_score($data[1]);
		$outfh->print(join(" ", @data) . "\n") if (defined $data[1] and $outfh);
		$count++;
		process_step_stats($data[1]) if $doStats;
	} 
	
	# a fixedStep data line
	elsif (scalar @data == 1) { 
		$data[0] = process_score($data[0]);
		$outfh->print($data[0] . "\n") if (defined $data[0] and $outfh);
		$count++;
		process_step_stats($data[0]) if $doStats;
	}
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
	my $mean  = sprintf("%.05f", $stats->{sumData} / $stats->{count});
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



# process
sub process_score {
	my $v = shift; # score
	if ($doNull and $v =~ /^(?:n.?[na])|(?:\-?inf)/i) {$v = 0}
	if ($deLogValue) {$v = $deLogValue ** $v}
	if ($doFlip) {$v *= -1}
	if ($multiplyValue) {$v *= $multiplyValue}
	if ($addValue) {$v += $addValue}
	if ($logValue) {$v = $v == 0 ? 0 : log($v) / $logValue}
	if ($doMin and $v < $minValue) {$v = $minValue}
	if ($doMax and $v > $maxValue) {$v = $maxValue}
	if ($places) {$v = sprintf($places, $v)};
	return undef if ($noZeroes and $v == 0);
	return $v;
}

sub process_bdg_stats {
	return unless defined $_[3];
	for ($_[1] + 1 .. $_[2]) {
		$stats->{count} += 1;
		$stats->{sumData} += $_[3];
		$stats->{sumSquares} += $_[3] ** 2;
		$stats->{minVal} = $_[0] if not defined $stats->{minVal};
		$stats->{maxVal} = $_[0] if not defined $stats->{maxVal};
		$stats->{minVal} = $_[3] if $_[3] < $stats->{minVal};
		$stats->{maxVal} = $_[3] if $_[3] > $stats->{maxVal};
	}
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



