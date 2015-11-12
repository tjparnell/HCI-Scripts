#!/usr/bin/env perl

use strict;
use Bio::ToolBox::Data;

unless (@ARGV) {
	die "Usage: $0 <file1.annot.txt> <file2.annot.txt> ...\n";
}

foreach my $file (@ARGV) {

	my $Data = Bio::ToolBox::Data->new(file => $file) or 
		die "unable to open $file!\n";
	
	my $ni = $Data->find_column('NORMAL');
	my $ti = $Data->find_column('TUMOR');
	my $nfi = $Data->add_column('NORMAL_AltFrequency');
	my $tfi = $Data->add_column('TUMOR_AltFrequency');
# 	$Data->reorder_column( (0..59), 61, 60, 62);
	
	$Data->iterate( sub {
		my $row = shift;
		
		# normal frequency
		if ($row->value($ni) =~ /:(\d+),(\d+)$/) {
			$row->value($nfi, sprintf("%.1f%%", ($2 / ($1 + $2)) * 100 ) );
		}
		elsif ($row->value($ni) =~ /:(\d+),(\d+),(\d+)$/) {
			$row->value(
				$nfi, sprintf("%.1f%%,%.1f%%", 
				( ($2 / ($1 + $2 + $3)) * 100 ), 
				( ($3 / ($1 + $2 + $3)) * 100 ) ) 
			);
		}
		else {
			warn sprintf("problem with line %s:%s\n", $row->seq_id, $row->start);
			$row->value($nfi, 0);
		}

		# tumor frequency
		if ($row->value($ti) =~ /:(\d+),(\d+)$/) {
			$row->value($tfi, sprintf("%.1f%%", ($2 / ($1 + $2)) * 100 ) );
		}
		elsif ($row->value($ti) =~ /:(\d+),(\d+),(\d+)$/) {
			$row->value(
				$tfi, sprintf("%.1f%%,%.1f%%", 
				( ($2 / ($1 + $2 + $3)) * 100 ), 
				( ($3 / ($1 + $2 + $3)) * 100 ) ) 
			);
		}
		else {
			warn sprintf("problem with line %s:%s\n", $row->seq_id, $row->start);
			$row->value($tfi, 0);
		}
	} );
	
	$Data->save;
	print "updated $file\n";
}

