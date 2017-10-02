#!/usr/bin/perl -w

use strict;
use warnings;
use POSIX ":sys_wait_h";
use Getopt::Long;
use Cwd;
use Bio::SearchIO;
use Bio::SeqIO;


#my $cds_dir = "/storage/RIBO_runs/RIBO_Bacteria/ANALYSIS/Conservation/1.Data/Protein";
my $cds_dir = "/data/elvis/REPARATION/ANALYSIS/Conservation/1.Data/Protein";

my $out_cds = "1.others_pep.fa";

#open(OUT, ">$out_cds")  or die $!;
opendir(DIRECTORY, $cds_dir) or die $!;
my $exit = 0;
while (my $file = readdir(DIRECTORY)) { 

	if($file=~/^\./){next;}
	my $file_in = path($file, $cds_dir);
	my ($header) = split(/\.([^\.]+)$/, $file);
	#$header =~	s/\.cds\.all//;
	$header =~	s/\.pep\.all//;
	print "$header\n";
	#my $cds_seq = read_fasta($file_in);

	#foreach my $id (keys %$cds_seq) {
		#print OUT ">".$header."#".$id."\n".uc($cds_seq->{$id})."\n";
	#}

	#last if ($exit > 1); $exit++;
}
close DIRECTORY;
#close OUT;

#system("usearch -cluster_fast z_bacteria_protein_oma_based.fa -id 0.9 -centroids z_bacteria_protein_oma_based_centroids.fa -uc clusters.uc");



##########################
#	Subs
##########################



sub read_fasta {

	my $file = $_[0];
	my $cdna = {};

	my $in  = Bio::SeqIO->new(-file => $file, -format => "fasta");
	while(my $seqs = $in->next_seq) {
		my $id  = $seqs->display_id;	
		my $seq = $seqs->seq;
		my $desc = $seqs->desc;

		#my ($chr) = $desc =~ /cds:annotated."?([^"]+)"gene/;
		my ($desc_info) = $desc =~ /pep:annotated (.*) gene:/;

		next unless($desc_info);

		my @info = split ':', $desc_info;
		next unless ($info[5]);
		my $strand = ($info[5] =~ '-') ? '-': '+';
		#print "$desc_info\t$strand\n"; exit;


		$id = $id."#".$info[2]."#".$info[3]."#".$info[4]."#".$strand;
		$cdna->{$id} = $seq;
	}
	return $cdna;
}


	# Check required parameters
sub uninitialized_param {
	
	my ($v) = @_;
	not ( defined($v) and length $v );
}

	# Handle path to files
sub path {

	my ($file, $dir) = @_;
	return (File::Spec->catfile( $dir, $file));
}
