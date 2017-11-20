#!/usr/bin/perl -w

####################################################################################################
#
#		Sarah B. Kingan
#		Pacific Biosciences
#		28 March 2017
#
#		Title: HiC_chaff_fasta.pl
#
#		Project: Aedes
#	
#		Input: 	1. tiling file from Olga
#				2. Unzip fragmented assembly
#				3. Unzip fragmented index file
#
#		Output: fasta file with sequences not included in scaffolded assembly
#			
#
####################################################################################################

use strict;

my $usage = "HiC_chaff.pl tiling.txt asm.fa asm.fa.fai\n";


# tiling file
my $tiling_file = shift(@ARGV) or die $usage;

# assembly file
my $asm_file = shift(@ARGV) or die $usage;

# index file
my $index_file = shift(@ARGV) or die $usage;

# contig length hash
my %contig_length_hash;
open (INDEX, $index_file);
while (my $line = <INDEX>) {
	my @line_array = split("\t", $line);
	$contig_length_hash{$line_array[0]} = $line_array[1];
}

# zero out fasta file
system('> tmp.fa');

# define variables
my $contigID = ''; # original ID
my $contigID1 = ''; # reformated for samtools faidx
my @ID_array = (); # to remove primary and secondary
my $ori = 0;		
my $rightClipLength = 0; # length clipped
my $leftClipLength = 0; # length clipped
my $start = 0;
my $end = 0;

open(TILING, "<", $tiling_file) or die("Can't open file");
my @lines = <TILING>;
close(TILING);
for (my $i = 1; $i<(scalar(@lines)-1); $i++) {
	my $line = $lines[$i];
	my @line_array = split("\t", $line);
	my $next_line = $lines[$i+1];
	my @next_line_array = split("\t", $next_line);
	unless ($line =~ /contig/) {
		@ID_array = split(" ", $line_array[1]);
		$contigID = $ID_array[0];
		$contigID1 = $contigID;
        $contigID1 =~ s/\|/\\\|/g; # format for samtools
        $ori = $line_array[0];
	}
	if (scalar(@line_array) == 2) { # first contig in mega contig
		if ($next_line =~ /contig/) { 
			# single contig megacontig-nothing clipped		
		}
		else { # first contig in multi mega contig
			$rightClipLength = $next_line_array[4] - $next_line_array[3];
			if ($rightClipLength > 200) {
				if ($ori == 0) {
					$start = $contig_length_hash{$contigID} - $rightClipLength + 1;
					$end = $contig_length_hash{$contigID};
				}
				elsif ($ori == 1) {
					$start = 1;
					$end = $start + $rightClipLength - 1;
				}
				system("samtools faidx $asm_file $contigID1:$start-$end >> tmp.fa");	
			}
			else {
				print $contigID, "\t", $rightClipLength, "\n";
			}
		}	
	}
	if (scalar(@line_array) == 8) {
		# left clip
		$leftClipLength = $line_array[6];
		if ($leftClipLength > 200) {
			if ($ori == 0) {
				$start = 1;
				$end = $leftClipLength;
			}
			elsif ($ori == 1) {
				$end = $contig_length_hash{$contigID};
				$start = $contig_length_hash{$contigID} - $leftClipLength + 1;
			}	
			system("samtools faidx $asm_file $contigID1:$start-$end >> tmp.fa");	
		}
		else {
			print $contigID, "\t", $leftClipLength, "\n";
		}
		
		# right clip
		if ($next_line =~ /contig/) { 
			# last contig in megacontig, no right clip
		}
		else {
			$rightClipLength = $next_line_array[4] - $next_line_array[3];
			if ($rightClipLength > 200) {
				if ($ori == 0) {
					$start = $contig_length_hash{$contigID} - $rightClipLength + 1;
             		$end = $contig_length_hash{$contigID};
				}
				elsif ($ori == 1) {
					$start = 1;
               		$end = $start + $rightClipLength - 1;
				}
				system("samtools faidx $asm_file $contigID1:$start-$end >> tmp.fa");	
			}
			else {
				print $contigID, "\t", $rightClipLength, "\n";
			}
		}
	}
}
			
