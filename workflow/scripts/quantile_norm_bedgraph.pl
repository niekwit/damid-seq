#!/usr/bin/perl -w
use strict;
use warnings;
use diagnostics;

$|++;

print STDERR "\nQuantile normalisation [filenames]\n\n";

# Load variables from the commandline
@ARGV || die "Need at least one .bedgraph file to process.\n\n";

# Globals
my @files;
my @probes;
my @scores;
my @scores_sorted;
my %data;

# Load files
my $counter=1;
foreach (@ARGV) {
	load_gff($_, $counter);
	$counter++;
}

my $num_files = @files;

# Re-format data
data_reformat();

# Normalisation
quant_norm();

# Write output
print STDERR "Writing output files ...              \n";
my $fn=1;
foreach my $f (@files) {
	write_output($f, $fn);
	$fn++;
}

print STDERR "All done.\n\n";

sub load_gff {
	my ($fn, $counter) = @_;
	push @files, $fn;															##Currently loaded file is stored in @file-array.
	print STDERR "Reading input file: $fn ...\n";
	open (IN, "<$fn") || die "Unable to read $fn: $!\n";

	my $count;
	foreach (<IN>) {
		$count++;
		my ($chr, $start, $end, $score) = split('\t');
		next unless $start;
		chomp $score;
		$chr=~s/chr//g;

		print STDERR "Read $count lines ...\r" if ($count%10000 == 0);

		# Final method:  setup hash keyed by chr and start pos, then store end pos at [0]
		# All scores follow.  Allows for uneven arrays, mismatching probes, the works.


		##20171127_script_changed_due_to_deprecated_defined()-command
		##-----------------------------------------------------------
		if ($counter == 1){
			push @{$data{$chr}{$start}},$end;
			push @{$data{$chr}{$start}},$score;
		} else {
			push @{$data{$chr}{$start}},$score;
		}

		##20171127_unless-control-structure_inactivated
		##---------------------------------------------
		#unless ( defined(@{$data{$chr}{$start}}) ){							##First time the probe is identified (probe=combination of chr,start,end) end will be stored.
		#	push @{$data{$chr}{$start}},$end;
		#}
		#push @{$data{$chr}{$start}},$score;									##For every '*.gff'-file, the score for the corresponding probe is stored at the end of the array.
	}
	close (IN);
}

sub data_reformat {
	my $p=0;
	foreach my $chr (sort keys %data) {											##Unroll every chromosome separately.
		foreach my $start (sort {$a <=> $b} keys $data{$chr}) {					##Loop through the start-positions of all lines and sort them.
			# exclude probes not covered in all data sets
			unless (@{ $data{$chr}{$start}} > $num_files) {						##If less scores are detected than the amount of files from the CLI, the probe was not uniformly identified throughout all files.
				# > num_files as first array datum is the end pos!
				print STDERR "Probe excluded: $chr:$start\n";
				next;
			}

			my @tempdata = @{ $data{$chr}{$start}};
			my $end = shift @tempdata;											##Separate the end_position from the beginning of the array, leaving only the scores.

			$probes[$p] = [$chr, $start, $end];									##Build an array of an array.
			@{ $scores[$p]} = ($p,@tempdata);									##Save scores together with number as identifier.
			#~ print "@{$scores[$p]}\n";
			$p++;
		}
	}
}

sub quant_norm {
	# Sorting
	my $total = $#scores;
	print STDERR "Sorting arrays ...       \n";
	foreach my $f (0 .. $#files) {												##Loop through the entire set of files as stored in the @file-array.
		my $g = $f+1;															##Add one in order for the first file to be numbered as '1' and not '0'.

		# sort the array on scores
		print STDERR "[Processing file $g] Sorting array ...                                  \r";
		@scores = sort { $a->[$g] <=> $b->[$g] } @scores;						##Order of scores per probe, i.e., @{ $scores[$p]}, is in line with the order of files in @file.
																				##All scores throughout one file will be resorted while still keeping the same position in the scores-array for every probe.
		# add each score to @scores_sorted (more cumbersome this way, takes more time, but greatly saves on memory)
		print "[Processing file $g] Adding elements to sorted scores ... \r";
		foreach my $r (0 .. $#scores) {
			if ($r%10000 == 0) {
				my $pc = sprintf("%0.2f", $r/$total*100);
				print STDERR "[Processing file $g] Adding elements to sorted scores: $pc% complete ...\r";
			}
			$scores_sorted[$r][$f] = $scores[$r][$g];							##Sorted array-of-arrays is stored by probe- and file-coordinates.
		}
	}

	# calculate rank scores
	print STDERR "Calculating and assigning rank scores ...                                \n";
	foreach my $r (0 .. $#probes) {												##Loop vertically over all the probes/GATC-fragments.
		if ($r%10000 == 0) {
			my $pc = sprintf("%0.2f", $r/$total*100);
			print STDERR "Processing: $pc% complete ...                        \r";
		}

		# Calculate rank averages
		my $rank_sum;
		foreach my $f (0 .. $#files) {											##Loop horizontally over all the sorted scores from the various files.
			$rank_sum += $scores_sorted[$r][$f];								##Sum up all the scores for that particular probe/GATC-fragment.
		}
		my $rank_score = $rank_sum/$num_files;									##Form the average over all scores in one probe.

		$scores_sorted[$r][$num_files] = $rank_score;							##Enhance array-of-arrays by another layer. Store rank-sum according to probe-coordinates.
	}

	# replace scores in original array
	print STDERR "Replacing scores in original arrays                        \n";
	foreach my $f (0 .. $#files) {
		my $g = $f+1;

		# sort the array on scores
		print STDERR "[Processing file $g] Sorting array ...                                  \r";
		@scores = sort { $a->[$g] <=> $b->[$g] } @scores;

		foreach my $r (0 .. $#scores) {
			if ($r%10000 == 0) {
				my $pc = sprintf("%0.2f", $r/$total*100);
				print STDERR "[Processing file $g] Replacing scores: $pc% complete ...                 \r";
			}
			$scores[$r][$g]=$scores_sorted[$r][$num_files];						##Replace the individual original scores distributed in probe-file-ccordinates with the rank-score.
		}
	}

	# Resort scores back to probes
	print STDERR "Re-sorting for output ...                                     \n";
	@scores = sort { $a->[0] <=> $b->[0] } @scores;								##Sort array back according to original probe-order.
}

sub write_output {
	my ($f, $fn) = @_;
	my $f_out = $f.".quant.norm.bedgraph";
	open (OUT, ">$f_out") || die "Cannot open $f_out for writing: $!\n";
	foreach my $p (0 .. $#probes) {
		my ($chr, $start, $end) = @{$probes[$p]};
		my $score = $scores[$p][$fn];
		$score = sprintf("%0.3f", $score);
		print OUT join("\t", $chr, $start, $end, $score), "\n";
	}
}

print "/nAll doneA.../n";
