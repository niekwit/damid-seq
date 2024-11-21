#!/usr/bin/perl -w

# Code adapted from:
# https://github.com/owenjm/damidseq_pipeline/blob/master/damidseq_pipeline
# Lines 786-911

use strict;
use warnings;
use diagnostics;

my @gatc_simple;
my %gatc;
my %gatc_index;

#my ($bam_in, $bam_out, %gatc_file, $lenvalue, $qvalue) = @ARGV;
my $bam_in = $ARGV[0];
my $bam_out = $ARGV[1];
my $gatc_file = $ARGV[2];
my $lenvalue = $ARGV[3];
my $qvalue = $ARGV[4];
my $logfile = $ARGV[5];

my %vars = (
    'bam' => $bam_in,
    'bam_out' => $bam_out,
    'gatc_frag_file' => $gatc_file,
    'len' => $lenvalue, 
    'q' => $qvalue,
    'samtools_path' => 'samtools ',
    'keep_original_bams' => 1,
    'log' => $logfile
);


sub load_gatc_frags {
	printout("\n*** Reading GATC file ***\n");
	
	if ($vars{'gatc_frag_file'} =~ m/\.gz$/) {
		# gzipped file
		open (GATC, '-|', "gunzip -c $vars{'gatc_frag_file'}") || die "Error: cannot open GATC file $vars{'gatc_frag_file'}: $!\n\n";
	} else {
		open (GATC, '<', $vars{'gatc_frag_file'}) || die "Error: cannot open GATC file $vars{'gatc_frag_file'}: $!\n\n";
	}
		
	foreach (<GATC>) {
		my ($chr, $source, $type, $start, $end, $score, $b, $c, $name) = split('\t');
		my $mid = ($start+$end)/2;
		push (@gatc_simple, [$chr, $mid]);
		push @{$gatc{$chr}}, $mid;
	}
	close GATC;
	
	print STDERR "  Sorting ...\n\n";
	my %tmp;
	foreach my $k (keys %gatc) {
		@{$tmp{$k}} = sort { $a <=> $b } @{$gatc{$k}};
	}
	%gatc = %tmp;
	
	# Build gatc_index lookup hash for new read extension method
	my $bin = 200;
	foreach my $k (keys %gatc) {
		
		my @tmp = @{ $gatc{$k} };
		
		foreach my $i (0 .. $#tmp) {
			my $site = $tmp[$i];
			my $bin = int($site/200)+1;
			$gatc_index{$k}{$bin} ||= $i;
		}
		
		my $last = 1;
		foreach my $i (1 .. max(keys %{ $gatc_index{$k}})) {
			$gatc_index{$k}{$i}||=$last;
			$last = $gatc_index{$k}{$i};
		}
	}
}


sub extend_reads_gatc {
    printout("*** Extending reads up to $vars{'len'} bases ***\n");

    printout("Reading input file: $vars{'bam'} ...\n");
    open (IN, "samtools view -h $vars{'bam'} |") || die "Unable to read $vars{'bam'}: $!\n";
    
    printout("  Processing data ...\n");
    open (OUT, "| samtools view -Shb | samtools sort -o $vars{'bam_out'}");
    
    my $c=0;
    my $seqs;
    my @non_chrs;
    while (<IN>) {
        
        if (/^\@/) {
            print OUT;
            next;
        }
        
        chomp;
        
        my ($qname, $flag, $chr, $pos, $mapq, $cigar, $rnext, $pnext, $tlen, $seq, $qual) = split('\s+');
        
        $c++;
        unless ($c%100000) {
            my $lp = sprintf("%0.2f",$c / 1000000);
            print "  $lp"."M reads processed ...\r";
        }
        
        next unless $seq;
        next unless $mapq>=$vars{'q'};
        
        unless ($gatc{$chr}) {
            # We do the missing chromosome check here if we're extending to GATCs
            push @non_chrs, $chr;
            next;
        }
        
        $seqs++;
        
        my $readlen;
        # used to use [M|I|S|=|X], but this is length of seq, not length of template
        foreach my $cig ($cigar =~ m/(\d+)[M|D|N|=|X]/g) {
            $readlen += $cig;
        }
        
        $cigar="$vars{'len'}M";
        
        # find next GATC after reads
        
        # It takes far too long to grep an entire chromosome's worth of GATC fragments for GATC sites internal to the extension length.  The following is an inelegant but very effective solution for speeding the process up: essentially, GATC fragments are binned into a hash when read in, and the hash is used as a lookup table.  There may be issues with long extension lengths, but it's fairly versatile.
        if ($flag & 16) {
            
            ### minus strand
            my $fiveprime_search = $pos - ($vars{'len'} - $readlen);
            my $threeprime_search = $pos;
            
            my $fp_bin = int($fiveprime_search/200)+1;
            my $tp_bin = int($threeprime_search/200)+2;
            
            my $fp_index = $gatc_index{$chr}{$fp_bin} || $#{ $gatc{$chr}}-2;
            my $tp_index = $gatc_index{$chr}{$tp_bin} || $#{ $gatc{$chr}}-2;
            
            my @search_slice = @{ $gatc{$chr}}[$fp_index .. $tp_index];
            
            my @hits = grep {$_ > $fiveprime_search && $_ < $threeprime_search} @search_slice;
            
            if (@hits) {
                # an intervening GATC
                my $closest = (sort {$a <=> $b} @hits)[$#hits];
                my $revised_len = $pos+$readlen - $closest;
                $cigar = $revised_len."M";
                $pos = $closest;
            } else {
                # no Great GATC
                $pos -= $vars{'len'} - $readlen;
                $pos = max(1,$pos);
            }
        } else {
            ### plus strand
            my $fiveprime_search = $pos + $readlen;
            my $threeprime_search = $pos + $vars{'len'};
                                
            my $fp_bin = int($fiveprime_search/200)+1;
            my $tp_bin = int($threeprime_search/200)+2;
            
            my $fp_index = $gatc_index{$chr}{$fp_bin} || $#{ $gatc{$chr}}-2;
            my $tp_index = $gatc_index{$chr}{$tp_bin} || $#{ $gatc{$chr}}-2;
            
            my @search_slice = @{ $gatc{$chr}}[$fp_index .. $tp_index];
            
            my @hits = grep {$_ > $fiveprime_search && $_ < $threeprime_search} @search_slice;
            
            if (@hits) {
                # an intervening GATC
                my $closest = (sort {$a <=> $b} @hits)[0];
                my $revised_len = $closest - $pos;
                $cigar = $revised_len."M";
            } else {
                # no Great GATC
                $pos = max(1,$pos);
            }
        }
        
        $seq = "*"; # We're extending reads and we have no sequence information for the extension ...
        $qual= "*"; #  ... ditto for quality.  Thankfully the SAMfile spec allows for no sequence or quality information.
        
        print OUT join("\t", $qname, $flag, $chr, $pos, $mapq, $cigar, $rnext, $pnext, $tlen, $seq, $qual), "\n";
        
    }
    close IN;
    close OUT;
    
    if (@non_chrs) {
        my $missing = join("\n    ", sort(uniq(@non_chrs)));
        print STDERR "  Warning: alignment contains chromosome identities not found in GATC file (see log file for details; this is normal if unmapped assemblies and the mitochondrial genome have been excluded from the GATC file.)\n";
        print LOG "  Warning: alignment contains chromosome identities not found in GATC file:\n    $missing\n\n";
    }
    
    printout("  Seqs extended (>q".$vars{'q'}.") = $seqs\n\n");
    #unlink($vars{'bam'}) unless $vars{'keep_original_bams'};
}

sub printout {
	my $s = shift;
	print STDERR $s;
	print LOG $s;
}

sub init_log_file {
    #my $log = $vars{'log'};
    open (LOG, '>', $vars{'log'});
    print LOG "Command-line values: @ARGV\n\n";
}

sub max {
    my ($max, @vars) = @_;
	my $index=0;
	$max||=0;
    for my $i (0..$#vars) {
        ($max, $index) = ($vars[$i], $i+1) if $vars[$i] > $max;
    }
    return ($index, $max);
}

# Log file
init_log_file();

# Load GATC fragments data
load_gatc_frags();

# Extend reads
extend_reads_gatc();
