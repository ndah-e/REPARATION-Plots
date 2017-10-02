#!/usr/bin/perl -w
	
use strict;
use warnings;
use Bio::SeqIO;

use Getopt::Long;
use List::Util qw/shuffle/;



my ($genome, $file, $fasta_file, $prefix) = @ARGV;

my %translationHash = 
	(GCA => "A", GCG => "A", GCT => "A", GCC => "A", GCN => "A",
     TGC => "C", TGT => "C",
     GAT => "D", GAC => "D",
     GAA => "E", GAG => "E",
     TTT => "F", TTC => "F",
     GGA => "G", GGG => "G", GGC => "G", GGT => "G", GGN => "G",
     CAT => "H", CAC => "H",
     ATA => "I", ATT => "I", ATC => "I",
     AAA => "K", AAG => "K",
     CTA => "L", CTG => "L", CTT => "L", CTC => "L", CTN => "L", TTA => "L", TTG => "L",
     ATG => "M",
     AAT => "N", AAC => "N",
     CCA => "P", CCT => "P", CCG => "P", CCC => "P", CCN => "P",
     CAA => "Q", CAG => "Q",
     CGA => "R", CGG => "R", CGC => "R", CGT => "R", CGN => "R",
     AGA => "R", AGG => "R",
     TCA => "S", TCG => "S", TCC => "S", TCT => "S", TCN => "S",
     AGC => "S", AGT => "S",
     ACA => "T", ACG => "T", ACC => "T", ACT => "T", ACN => "T",
     GTA => "V", GTG => "V", GTC => "V", GTT => "V", GTN => "V",
     TGG => "W",
     TAT => "Y", TAC => "Y",
     TAG => "*", TAA => "*", TGA => "*");

my @AA = ("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y");


# read genome and RPF occupancy files
my $genomes = read_fasta($genome);

my $ORFs = read_file($file);

#calculate_codon_freq($ORFs);
calculate_codon_freq_2nd($ORFs,$fasta_file);

#calculate_dicodon_freq($ORFs);
#system("Rscript /data/elvis/REPARATION/github/ANALYSIS/scripts/plots_codon_usage.R $file $merge_file");



##################
##	SUBS
##################


sub calculate_codon_freq_2nd {

	my $ORFs = $_[0];
	my $fasta_file = $_[1];

	my $annotated = {};
	my $novel = {};

	my $count_anno = 0;
	my $count_novel = 0;
	my $count_random = 0;


	foreach my $ORF (keys %$ORFs) {

		my $start = $ORFs->{$ORF}->{start};
		my $stop = $ORFs->{$ORF}->{stop};
		my $region = $ORFs->{$ORF}->{region};
		my $strand = $ORFs->{$ORF}->{strand};
		my $ORF_type = $ORFs->{$ORF}->{ORF_type};
		my $length = $ORFs->{$ORF}->{len};
		my $nterm = $ORFs->{$ORF}->{nterm};

		#next if ($ORF_type eq 'Extension');
		#next if ($ORF_type eq 'Truncation');

		my $dna_seq = ($strand eq '+') ? substr($genomes->{$region}, $start - 1,$length) : revdnacomp(substr($genomes->{$region}, $start-1,$length));
		for (my $i = 3; $i < 6; $i = $i + 3) {
			last if (length(substr($dna_seq, $i, 3)) < 3);
			my $aa = $translationHash{uc(substr($dna_seq, $i, 3))};
			if ($aa ne '*') {
				if ($ORF_type eq 'Annotated') {
					#print "$ORF\t$aa\n";
					$annotated->{$aa}++;
					$count_anno++;
				} else {
					$novel->{$aa}++;
					$count_novel++;
				}
			}
		}
	}

	my $random = average_aa($fasta_file);

	my $file = $prefix.".txt";
	open(F, ">".$file);
	print F"aa\tAnnotated\tNovel\tRandom\n";
	for my $aa (@AA) {
		my $anno = ($annotated->{$aa}) ? $annotated->{$aa}/$count_anno: 0;
		my $novel_c = ($novel->{$aa}) ? $novel->{$aa}/$count_novel: 0;
		#my $random_c = ($random->{$aa}) ? $random->{$aa}/$count_random: 0;
		my $random_c = ($random->{$aa}) ? $random->{$aa}: 0;

		print F "$aa\t$anno\t$novel_c\t$random_c\n";
	}
	close F;

	system("Rscript /data/elvis/REPARATION/ANALYSIS/Di_codon/dicodon.R $file $prefix");
}


sub average_aa {

	my $file = $_[0];

	my $random = {};
	my $fasta = read_fasta($file);

	my $count = 0;
	foreach my $id (keys %$fasta) {
		my @aa= split('', $fasta->{$id});
		foreach my $a (@aa) {
			$count++;
			$random->{$a}++;
		}
	}

	foreach my $a (keys %$random) {
		$random->{$a} = $random->{$a}/$count;
	}

	return $random;
}


sub read_file {

	my $file = $_[0];
	my $ORFs = {};

	open (F, $file) or die  "Error reading file: $file";
	while (<F>) {
		chomp $_;
		if (/^orf_id/) {next;}
		if (/^ORF_locus/) {next;}
		my @line 			= (split '\t', $_);

		my $ORF 	= $line[0];
		my $strand 	= $line[1];
		my $length 	= $line[2];

		my $start_codon = $line[4];
		my $ref = $line[5];
		my $nterm = $line[6];
		my $dist = $line[7];

		my $ORF_type 	= $line[3];
		#my $ORF_type 	= $line[10];

		if ($ORF_type eq 'Extension' and abs($dist) < 9) {next}
		#if ($ORF_type eq 'Truncation' and abs($dist) < 9) {next}

		my ($start) = $ORF =~ /:(\d+)-/;
		my ($stop) = $ORF =~ /-(\d+)$/;
		my ($region) = $ORF =~ /^(.*):/;

		my $gene = ($strand eq "+") ? $region.":FAM+".$stop: $region.":FAM-".$start;

		$ORFs->{$ORF}->{start}	= $start;
		$ORFs->{$ORF}->{stop}	= $stop;
		$ORFs->{$ORF}->{region} = $region;
		$ORFs->{$ORF}->{strand} = $strand;
		$ORFs->{$ORF}->{ORF_type}= $ORF_type;
		$ORFs->{$ORF}->{len} 	= $length;
		$ORFs->{$ORF}->{nterm} 	= $nterm;
		$ORFs->{$ORF}->{dist} 	= $dist;
	}
	close F; 

	return $ORFs;
}


sub generate_seq {

	my $length = $_[0];

	$length = $length/3;

	rand($length);
	#rand($length);

	#my @chars = ("C","G","A","T");
	my @chars = @AA;
	my $string;
	$string .= $chars[rand @chars] for 1..$length;

	#print "$string\n\n";
	return $string;
}

sub translate {

	my $start = $_[0];
	my $stop = $_[1];
	my $strand = $_[2];
	my $region = $_[3];
	my $ORF = $_[4];

	my $length = $stop - $start + 1;
	my $dna_seq = ($strand eq '+') ? substr($genomes->{$region}, $start - 1,$length) : revdnacomp(substr($genomes->{$region}, $start-1,$length));

	my $aa_seq = "";
	for (my $i = 0; $i <= (length($dna_seq) - 3); $i = $i + 3) {
		my $codon = uc(substr($dna_seq, $i, 3));
		my $aa = $translationHash{$codon};
		last if ($aa eq '*');
		$aa_seq = $aa_seq.$aa;
	}
	return $aa_seq;
}

sub read_fasta {
	my $file = $_[0];
	my $cdna = {};

	my $in  = Bio::SeqIO->new(-file => $file, -format => "fasta");
	while(my $seqs = $in->next_seq) {
		my $id  = $seqs->display_id;	
		my $seq = $seqs->seq;
		my $desc = $seqs->desc;
		$cdna->{$id} = uc($seq);
	}
	return $cdna;
}


sub revdnacomp {
    my $dna = shift;
    my $revcomp = reverse($dna);
    $revcomp =~ tr/ACGTacgt/TGCAtgca/;
    return $revcomp;
}




sub calculate_dicodon_freq {

	my $ORFs = $_[0];

	my $annotated = {};
	my $novel = {};
	my $random = {};

	my $count_anno = 0;
	my $count_novel = 0;
	my $count_random = 0;
	foreach my $ORF (keys %$ORFs) {

		my $start = $ORFs->{$ORF}->{start};
		my $stop = $ORFs->{$ORF}->{stop};
		my $region = $ORFs->{$ORF}->{region};
		my $strand = $ORFs->{$ORF}->{strand};
		my $ORF_type = $ORFs->{$ORF}->{ORF_type};
		my $length = $ORFs->{$ORF}->{len};

		next if ($ORF_type eq 'Extension');
		next if ($ORF_type eq 'Truncation');

		my $dna_seq = ($strand eq '+') ? substr($genomes->{$region}, $start - 1,$length) : revdnacomp(substr($genomes->{$region}, $start-1,$length));

		my $aa_seq = "";
		for (my $i = 0; $i <= (length($dna_seq) - 3); $i = $i + 3) {
			last if (length(substr($dna_seq, $i, 3)) < 3);
			my $aa = $translationHash{uc(substr($dna_seq, $i, 3))};
			if ($aa ne '*') {
				$aa_seq = $aa_seq.$aa;
			}
		}

		for (my $i = 0; $i <= (length($aa_seq) - 2); $i = $i + 2) {
			my $di_aa = substr($aa_seq,$i,2);
			if ($ORF_type eq 'Annotated') {
				$annotated->{$di_aa}++;
				$count_anno++;
			} else {
				$novel->{$di_aa}++;
				$count_novel++;
			}
		}

		# random sequences
		my $dna_seq_rev = ($strand eq '+') ? revdnacomp(substr($genomes->{$region}, $start-1,$length)) : substr($genomes->{$region}, $start - 1,$length);
		my $aa_seq_rev = "";
		for (my $i = 0; $i <= (length($dna_seq_rev) - 3); $i = $i + 3) {
			next if (length(uc(substr($dna_seq_rev, $i, 3))) < 3);
			my $aa = $translationHash{uc(substr($dna_seq_rev, $i, 3))};
			if ($aa ne '*') {
				$aa_seq_rev = $aa_seq_rev.$aa;
			}
		}

		for (my $i = 0; $i <= (length($aa_seq_rev) - 2); $i = $i + 2) {
			my $di_aa = substr($aa_seq_rev,$i,2);
			$random->{$di_aa}++;
			$count_random++;
		}
	}

	my $dicodons = {%$annotated,%$novel,%$random};
	my $file = $prefix.".txt";
	#open(F, ">".$file);
	print "dicodon\tAnnotated\tNovel\tRandom\n";
	for my $aa (keys %$dicodons) {
		my $anno = ($annotated->{$aa}) ? $annotated->{$aa}/$count_anno: 0;
		my $novel_c = ($novel->{$aa}) ? $novel->{$aa}/$count_novel: 0;
		my $random_c = ($random->{$aa}) ? $random->{$aa}/$count_random: 0;

		print "$aa\t$anno\t$novel_c\t$random_c\n";
	}
	close F;
}




sub calculate_codon_freq_old {

	my $ORFs = $_[0];

	my $annotated = {};
	my $novel = {};
	my $random = {};

	my $count_anno = 0;
	my $count_novel = 0;
	my $count_random = 0;


	foreach my $ORF (keys %$ORFs) {

		my $start = $ORFs->{$ORF}->{start};
		my $stop = $ORFs->{$ORF}->{stop};
		my $region = $ORFs->{$ORF}->{region};
		my $strand = $ORFs->{$ORF}->{strand};
		my $ORF_type = $ORFs->{$ORF}->{ORF_type};
		my $length = $ORFs->{$ORF}->{len};
		my $nterm = $ORFs->{$ORF}->{nterm};

		next if ($ORF_type eq 'Extension');
		next if ($ORF_type eq 'Truncation');

		my $dna_seq = ($strand eq '+') ? substr($genomes->{$region}, $start - 1,$length) : revdnacomp(substr($genomes->{$region}, $start-1,$length));
		for (my $i = 0; $i <= (length($dna_seq) - 3); $i = $i + 3) {
			last if (length(substr($dna_seq, $i, 3)) < 3);
			my $aa = $translationHash{uc(substr($dna_seq, $i, 3))};
			if ($aa ne '*') {
				if ($ORF_type eq 'Annotated') {
					print "$ORF\t$aa\n";
					$annotated->{$aa}++;
					$count_anno++;
				} else {
					#if ($nterm eq 'Y' or $nterm eq 'V') {
					#if ($nterm eq 'N') {
					$novel->{$aa}++;
					$count_novel++;
					#}
				}
			}
		}

		next if ($ORF_type ne 'Annotated');
		# random sequences = Reverse complement
		my $dna_seq_rev=  join('',shuffle(split(//,$dna_seq)));
		#my $dna_seq_rev = ($strand eq '+') ? revdnacomp(substr($genomes->{$region}, $start-1,$length)) : substr($genomes->{$region}, $start - 1,$length);
		my $aa_seq_rev = "";
		for (my $i = 0; $i <= (length($dna_seq_rev) - 3); $i = $i + 3) {
			next if (length(uc(substr($dna_seq_rev, $i, 3))) < 3);
			my $aa = $translationHash{uc(substr($dna_seq_rev, $i, 3))};
			if ($aa ne '*') {
				$random->{$aa}++;
				$count_random++;
			}
		}

		# random sequences = frame offset
		#my $dna_seq_rev = ($strand eq '+') ? substr($genomes->{$region}, $start - 1,$length) : revdnacomp(substr($genomes->{$region}, $start-1,$length));
		#for (my $i = 2; $i <= (length($dna_seq_rev) - 3); $i = $i + 3) {
		#	next if (length(uc(substr($dna_seq_rev, $i, 3))) < 3);
		#	my $aa = $translationHash{uc(substr($dna_seq_rev, $i, 3))};
		#	if ($aa ne '*') {
		#		$random->{$aa}++;
		#		$count_random++;
		#	}
		#}

		#my @rand_aa= split('', generate_seq($length));
		#for my $aa (@rand_aa) {
		#	$random->{$aa}++;
		#	$count_random++;
		#}
	}

	my $file = $prefix.".txt";
	open(F, ">".$file);
	print F"aa\tAnnotated\tNovel\tRandom\n";
	for my $aa (@AA) {
		my $anno = ($annotated->{$aa}) ? $annotated->{$aa}/$count_anno: 0;
		my $novel_c = ($novel->{$aa}) ? $novel->{$aa}/$count_novel: 0;
		my $random_c = ($random->{$aa}) ? $random->{$aa}/$count_random: 0;

		print F "$aa\t$anno\t$novel_c\t$random_c\n";
	}
	close F;

	system("Rscript /data/elvis/REPARATION/ANALYSIS/Di_codon/dicodon.R $file $prefix");
}

