#!/usr/bin/perl

# a script to run Java TreeView to create images

use strict;
use File::Spec;
use Getopt::Long;


# help
unless (@ARGV) {
print <<HELP;
  A wrapper script for printing TreeView images from cluster data files
  Usage:
    print_clusters.pl [options] file1.cdt file2.cdt ...
  Options:
    -c contrast, default 3
    -s scale width X height in pixels, default 2x1
    -f picture format, default png
    -j java path, default /usr/bin/java
    -t TreeView jar path, default ~/bin/TreeView-1.1.6r2-bin/TreeView.jar
    -m java memory, default 2G
  See online help at http://jtreeview.sourceforge.net/docs/JTVUserManual/ch02.html#ch2_flags
  for more information
HELP
exit;
}

# default options
my $contrast = 3;
my $scale = '2x1';
my $format = 'png';
my $java = '/usr/bin/java';
my $treeview = "$ENV{HOME}/bin/TreeView-1.1.6r2-bin/TreeView.jar";
my $memory = '2G';

# command line options
GetOptions(
	'c=s'   => \$contrast,
	's=s'   => \$scale,
	'f=s'   => \$format,
	'j=s'   => \$java,
	't=s'   => \$treeview,
	'm=s'   => \$memory,
) or die "bad options!\n";

foreach (@ARGV) {
	
	# files
	my $file = File::Spec->rel2abs($_);
	print "##### Generating image for $file\n";
	my $out = $file;
	$out =~ s/\.cdt$/.$format/;
	
	# command
	# java -jar ~/bin/TreeView-1.1.6r2-bin/TreeView.jar -x Dendrogram -r ./TSS_Rsc1_with_Rsc2_data.cdt -- -o ./TSS_Rsc1_with_Rsc2_data.png -c 1.5 -s 4x0.2 -f png 
	system(
		$java,
		'-Xmx'. $memory,
		'-jar',
		$treeview,
		'-x',
		'Dendrogram',
		'-r',
		$file,
		'--',
		'-o',
		$out,
		'-f',
		$format,
		'-s',
		$scale,
		'-c',
		$contrast,
	) == 0 or warn "##### $file didn't work!! #######";
	
	# check
	if (-e $out) {
		print "##### $out file generated!\n";
	}
}
