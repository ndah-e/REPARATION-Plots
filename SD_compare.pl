#!/usr/bin/perl -w
	
use strict;
use warnings;
#use diagnostics;

use Getopt::Long;
use Bio::SeqIO;


#######################################


my ($infile,$genome);
my $expname = "context";
my $REGION = 20;


GetOptions(
	'i=s' =>\$infile,
	'g=s' =>\$genome,
	'e=s' =>\$expname
);

my %params = (
	i => $infile,
	g => $genome,
);


my $genomes = read_fasta($genome);

my $ORFs = read_file($infile);
sequence_context($ORFs);


##################
##	SUBS
##################



sub sequence_context {

	my $ORFs = $_[0];
	
	my $file = $expname."_context.txt";
	open(F, ">".$file) or die  "Error creating file: $file\n";
	print F "ORF_type\tnterm\tcontext\tcontext_ref\n";
	
	foreach my $ORF (keys %$ORFs) {

		my $region = $ORFs->{$ORF}->{region};
		my $strand = $ORFs->{$ORF}->{strand};
		my $ORF_type = $ORFs->{$ORF}->{ORF_type};
		my $reference = $ORFs->{$ORF}->{reference};
		my $nterm = $ORFs->{$ORF}->{nterm};
		my $start = $ORFs->{$ORF}->{start};
		my $stop = $ORFs->{$ORF}->{stop};
		
		if ($ORF_type eq 'Truncation' or $ORF_type eq 'Extension') {

			#next if (abs($ORFs->{$ORF}->{dist_aTIS}) < 20);
			my ($start_ref) = $reference =~ /:(\d+)-/;
			my ($stop_ref) 	= $reference =~ /-(\d+)$/;
			
			my $seq = $genomes->{$region};
			if ($strand eq '+') {
				my $context = substr($seq, $start - $REGION - 1, $REGION + 3);
				my $context_ref = substr($seq, $start_ref - $REGION - 1, $REGION + 3);
				print F "$ORF_type\t$nterm\t$context\t$context_ref\n";
				
			} else {
				my $context = revdnacomp(substr($seq, $stop - 3, $REGION+3));
				my $context_ref = revdnacomp(substr($seq, $stop_ref - 3, $REGION+3));
				print F "$ORF_type\t$nterm\t$context\t$context_ref\n";
			}
		}
	}
	close F;
}


sub read_file {

	my $file = $_[0];
	my $ORFs = {};

	open (F, $file) or die  "Error reading file: $file";
	while (<F>) {
		chomp $_;
		if (/^orf_id/) {next;}
		
		my @line 		= (split '\t', $_);
		my $ORF 		= $line[0];
		my $gene		= $line[1];
		my $strand 		= $line[2];
		my $length 		= $line[3];
		my $start_codon = $line[4];
		
		my $ORF_type	= $line[16];
		my $dist_aTIS	= $line[17];
		my $reference	= $line[18];
		my $nterm		= $line[19];
		
		($ORFs->{$ORF}->{start}) 	= $ORF =~ /:(\d+)-/;
		($ORFs->{$ORF}->{stop}) 	= $ORF =~ /-(\d+)$/;
		my ($region) = $ORF =~ /^(.*):/;

		$ORFs->{$ORF}->{region} 	= $region;
		$ORFs->{$ORF}->{strand} 	= $strand;
		$ORFs->{$ORF}->{len} 		= $length;
		$ORFs->{$ORF}->{start_cdn} 	= $start_codon;
		#$ORFs->{$ORF}->{SD_seq} 	= $SD_seq;
		$ORFs->{$ORF}->{ORF_type} 	= $ORF_type;
		$ORFs->{$ORF}->{reference} 	= $reference;
		$ORFs->{$ORF}->{dist_aTIS}	= $dist_aTIS;
		$ORFs->{$ORF}->{nterm}		= $nterm;

	}
	close F; 

	return $ORFs;
}


sub revdnacomp {
    my $dna = shift;
    my $revcomp = reverse($dna);
    $revcomp =~ tr/ACGTacgt/TGCAtgca/;
    return $revcomp;
}



sub read_fasta {
	my $file = $_[0];
	my $cdna = {};

	my $in  = Bio::SeqIO->new(-file => $file, -format => "fasta");
	while(my $seqs = $in->next_seq) {
		my $id  = $seqs->display_id;	
		my $seq = $seqs->seq;
		my $desc = $seqs->desc;
		$cdna->{$id} = $seq;
	}
	return $cdna;
}



sub read_file_old {

	my $file = $_[0];
	my $ORFs = {};

	open (F, $file) or die  "Error reading file: $file";
	while (<F>) {
		chomp $_;
		if (/^orf_id/) {
			$ORFs->{gene}->{ORF}->{head} = $_;
			next;
		}
		my @line 			= (split '\t', $_);
		my $ORF 			= $line[0];
		my $gene			= $line[1];
		my $strand 			= $line[2];
		my $length 			= $line[3];
		my $start_codon 	= $line[4];
		my $ribo_count		= $line[5];
		my $ribo_rpkm		= $line[6];
		my $ribo_coverage	= $line[7];
		my $start_coverage	= $line[8];
		my $start_rpkm 		= $line[9];
		my $stop_rpkm 		= $line[10];
		my $orf_score		= $line[11];
		my $SD_score 		= $line[12];
		my $SD_seq 			= $line[13];
		my $prob			= $line[14];
		my $pc				= $line[15];
		my $rna_count 		= $line[16];
		my $rna_rpkm		= $line[17];
		my $TE				= $line[18];
		my $mmass			= $line[19];
		my $ORF_type		= $line[20];
		my $dist_aTIS		= $line[21];
		my $reference		= $line[22];
		my $nterm			= $line[23];
		
		($ORFs->{$gene}->{$ORF}->{start}) 	= $ORF =~ /:(\d+)-/;
		($ORFs->{$gene}->{$ORF}->{stop}) 	= $ORF =~ /-(\d+)$/;
		my ($region) = $ORF =~ /^(.*):/;

		$ORFs->{$gene}->{$ORF}->{orf_info} 	= $_;
		$ORFs->{$gene}->{$ORF}->{region} 	= $region;
		$ORFs->{$gene}->{$ORF}->{strand} 	= $strand;
		$ORFs->{$gene}->{$ORF}->{len} 		= $length;
		$ORFs->{$gene}->{$ORF}->{start_cdn} = $start_codon;
		$ORFs->{$gene}->{$ORF}->{start_rpkm}= $start_rpkm;
		$ORFs->{$gene}->{$ORF}->{ribo_rpkm} = $ribo_rpkm;
		$ORFs->{$gene}->{$ORF}->{ribo_count}= $ribo_count;
		$ORFs->{$gene}->{$ORF}->{coverage} 	= $ribo_coverage;
		$ORFs->{$gene}->{$ORF}->{SD_score} 	= $SD_score;
		$ORFs->{$gene}->{$ORF}->{SD_seq} 	= $SD_seq;
		$ORFs->{$gene}->{$ORF}->{orf_score} = $orf_score;
		$ORFs->{$gene}->{$ORF}->{prob} 		= $prob;
		$ORFs->{$gene}->{$ORF}->{pc} 		= $pc;
		$ORFs->{$gene}->{$ORF}->{rna_count} = $rna_count;
		$ORFs->{$gene}->{$ORF}->{rna_rpkm} 	= $rna_rpkm;
		$ORFs->{$gene}->{$ORF}->{TE} 		= $TE;
		$ORFs->{$gene}->{$ORF}->{ORF_type} 	= $ORF_type;
		$ORFs->{$gene}->{$ORF}->{reference} = $reference;
		$ORFs->{$gene}->{$ORF}->{nterm}		= $nterm;

	}
	close F; 

	return $ORFs;
}

