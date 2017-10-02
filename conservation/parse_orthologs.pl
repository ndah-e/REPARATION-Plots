#!/usr/bin/perl -w

use strict;
use warnings;

use Getopt::Long;
use Bio::SeqIO;
use Cwd;

my $genome_dir = "/storage/RIBO_runs/RIBO_Bacteria/ANALYSIS/Conservation/1.Data/genomes";
my $ortho_file = "/storage/RIBO_runs/RIBO_Bacteria/ANALYSIS/Conservation/SALT/Annotated/fasta/Results_Aug18/OrthologousGroups.csv";
my $work_dir = "";
my $ref_genome;
my $muscle = "/data/elvis/REPARATION/ANALYSIS/Conservation/scripts_n_tools/muscle3.8.31_i86linux64";
my $outprefix;
my $predicted;

my $REGION = 30;
my $MINORF = 30;
my $MINORTH = 5;

GetOptions(
	'i=s' =>\$ortho_file,
	'g=s' =>\$genome_dir,
	'r=s' =>\$ref_genome,
	'o=s' =>\$outprefix,
	'l=i' =>\$MINORF,
	'c=i' =>\$MINORTH,
	'w=s' =>\$work_dir,
	'p=s' =>\$predicted	# to remove from script
);

my %params = (
	i => $ortho_file,
	r => $ref_genome,
	#g => $genome_dir

);

my @invalid = grep uninitialized_param($params{$_}), keys %params;	
die "Not properly initialized: @invalid\n" if @invalid;

#if ($outprefix) {$outprefix = $outprefix."_"} else {$outprefix = ""}

my ($workdir, $tmp_dir) = check_working_dir($work_dir);
#print "$workdir\n$tmp_dir\n";

#my $genomes = read_fasta($genome);

#my $ORFs = read_file($predicted);

my $orthologs = get_ortholog($ortho_file);

process_orthologs($orthologs);


print "Site Conservation ccalculations completed\n";


##########################
#	Subs
##########################


sub process_orthologs {

	my $orthologs = $_[0];
	
	my $rate4site = {};
	my $count_seq = 0;
	foreach my $orf (keys %$orthologs) {

		my $orf_path = path($orf, $tmp_dir);

		open(F, ">".$orf_path.".fa") or die "can not creat file $orf_path.fa\n";

		my ($region) = $orf =~ /^(.*):/;
		my ($start) = $orf =~ /:(\d+)-/;
		my ($stop) 	= $orf =~ /-(\d+)#/;
		my ($strand) = $orf =~ /#(.)$/;
		my $length = $stop - $start + 1;

		my $orf_id = (split '#', $orf)[0];

		#next if (abs($ORFs->{$orf_id}->{dist}) <= 9);

		my $ref_genomes = read_fasta($ref_genome);

		my $orf_dna_seq = ($strand eq '+') ? substr($ref_genomes->{$region}, $start - 1,$length) : revdnacomp(substr($ref_genomes->{$region}, $start-1,$length));
		my $consensus = ($strand eq '+') ? substr($ref_genomes->{$region}, $start - $REGION - 1,$REGION) : revdnacomp(substr($ref_genomes->{$region}, $stop,$REGION));
		$orthologs->{$orf}->{consen} = $consensus;

		print F ">$orf\n$orf_dna_seq\n";

		foreach my $ortholog (keys %{$orthologs->{$orf}}) {

			next if ($ortholog eq "consen");
			my @orth_info = split "#", $ortholog;

			my $g_file = $orth_info[0];
			my $region = $orth_info[2];
			my $start = $orth_info[3];
			my $stop = $orth_info[4];
			my $strand = $orth_info[5];

			$g_file = $g_file.".dna.toplevel.fa";
			$g_file = path($g_file, $genome_dir);
			my $genomes = read_fasta($g_file);

			my $length = $stop - $start + 1;
			my $dna_seq = ($strand eq '+') ? substr($genomes->{$region}, $start - 1,$length) : revdnacomp(substr($genomes->{$region}, $start-1,$length));

			my $consensus = ($strand eq '+') ? substr($genomes->{$region}, $start - 1 - $REGION,$REGION) : revdnacomp(substr($genomes->{$region}, $stop,$REGION));

			if ($strand eq '+') {
				if ($start - 1 - $REGION < 0) {
					my $lenth_tmp = $REGION + ($start - 1 - $REGION);
					$consensus = substr($genomes->{$region}, 0,$lenth_tmp);
				}
			} else {
				if ($stop + $REGION > length($genomes->{$region})) {
					my $lenth_tmp = $REGION + (length($genomes->{$region}) - ($stop + $REGION));
					$consensus = substr($genomes->{$region}, $stop,$lenth_tmp);
				}
			}

			if (length($consensus) < $REGION) {
				for (my $i = 0; $i < $REGION - length($consensus); $i++) {
					$consensus = "-".$consensus;
				}
			} 

			$orthologs->{$orf}->{$ortholog} = $consensus;

			print F ">$ortholog\n$dna_seq\n";
		}
		close F;

		# perform muscle
		my $muscle_out = $orf_path."_algn.fa";
		system($muscle." -in ".$orf_path.".fa -out $muscle_out -fasta -quiet");

		# attach UTR
		my $alignment = parse_fasta($muscle_out);
		my $muscle_align = read_fasta($muscle_out);

		# get start position after alignment
		my $start_seq = substr($orf_dna_seq,0,30);
		my $align_start = find_orf_start($start_seq,$muscle_align->{$orf});

		my $file_utr = $orf_path."_algn_utr.fa";
		open (FU, ">".$file_utr) or die "error creating $file_utr\n";
		foreach my $count (sort {$a <=> $b} keys %$alignment) {
			my $id = (keys %{$alignment->{$count}})[0];
			my $utr = ($id eq $orf) ? $orthologs->{$id}->{consen} :$orthologs->{$orf}->{$id};
			print FU ">$id\n".$utr.$alignment->{$count}->{$id}."\n";
		}
		close FU;

		# run rate4site
		my $r4s_out = $orf_path."_r4s_file";
		system("rate4site -s $file_utr -o $r4s_out -a $orf -Mn");

		# parse rate4site result
		open(FR, $r4s_out) or die "can not creat file $r4s_out\n";
		while (<FR>) {

			chomp $_;
			next if ($_ =~ /^#/);
			$_ =~ s/^\s+//;
			next if ($_ =~ /^\s*$/);

			my @line = split '\t|\s+', $_;
			my $pos = $line[0] - $align_start - $REGION;
			my $score = -1*$line[2];

			$rate4site->{$pos}->{score} += $score;
			$rate4site->{$pos}->{count} += 1;
		}
		close FR;

		$count_seq++;
	}

	$outprefix =~ s/\s+//;
	open(F, ">".$outprefix."_output_agregate.txt") or die ;
	my $pos = "posA".$count_seq;
	print F "$pos\taverage_score\n";
	foreach my $pos (sort {$a <=> $b} keys %$rate4site) {
		my $average = $rate4site->{$pos}->{score}/$rate4site->{$pos}->{count};
		print F "$pos\t$average\n";
	}
	close F;

	print "\n\n**************\n\nNumber of sequences in meta gene plot $count_seq\n\n**************\n\n";

}



sub find_orf_start {

	my $start_seq = $_[0];
	my $sequence = $_[1];

	my @start = split '', $start_seq;
	my $search_string;
	foreach my $codon (@start){
		if ($search_string) {
			$search_string = $search_string."(-*)".$codon;
		} else {
			$search_string = $codon."(-*)";
		}
	}

	#$search_string = "(".$search_string.")";
	#print "$search_string\n";
	my $pos;
	if ($sequence=~ /($search_string)/){
		#$pos = pos($sequence)-length($1);
		#$pos = $-[0] + 1;
		$pos = $-[0];
	}

	return $pos;
}


sub get_ortholog {

	my $file = $_[0];

	my $orthologs = {};
	open(F, $file) or die "cannot read file $file\n";
	while (<F>) {

		next if (/^\s+/);
		chomp $_;

		my @line = split '\t', $_;
		$line[1] =~ s/\s+//g;
		my @orfs = split ',', $line[1];
		foreach my $orf (@orfs) {

			$line[2] =~ s/\s+//g;
			my @ortho = split ',', $line[2];
			next if (scalar(@ortho) < $MINORTH);

			foreach my $ortholog (@ortho) {
				if ($ortholog =~ /:/) {print "homolog $ortholog\n"; next}		# skip any homolog/ortholog in input file
				$orthologs->{$orf}->{$ortholog} = 1;
			}
		}
	}
	close F;

	return $orthologs;
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


sub parse_fasta {

	my $file = $_[0];
	my $fasta = {};

	my $count = 1;
	my $in  = Bio::SeqIO->new(-file => $file, -format => "fasta");
	while(my $seqs = $in->next_seq) {
		my $id  = $seqs->display_id;	
		my $seq = $seqs->seq;
		$fasta->{$count}->{$id} = $seq;
		$count++;
	}
	return $fasta;
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

sub revdnacomp {
    my $dna = shift;
    my $revcomp = reverse($dna);
    $revcomp =~ tr/ACGTacgt/TGCAtgca/;
    return $revcomp;
}

sub check_working_dir {

	my $work_dir = $_[0];

	my @months = qw( Jan Feb Mar Apr May Jun Jul Aug Sep Oct Nov Dec );
	my @days = qw(Sun Mon Tue Wed Thu Fri Sat Sun);
	my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime();
	my $tmp_dir;

	# Check temp directory
	if ($work_dir) {
		#$tmp_dir = $work_dir."/".$months[$mon]."_".$days[$wday]."_".$sec;
		$tmp_dir = $work_dir."/temp_".$months[$mon]."_".$mday;
		if (!-d $work_dir) {
			system("mkdir -p $work_dir");
			system("mkdir -p $tmp_dir");	# recreate tmp directory
		} else {
			system("rm -rf $work_dir");
			system("rm -rf $tmp_dir");
			system("mkdir -p $work_dir");
			system("mkdir -p $tmp_dir");
		}
	} else {
		$work_dir = getcwd();
		#$tmp_dir = $work_dir."/".$months[$mon]."_".$days[$wday].$sec;
		$tmp_dir = $work_dir."/".$months[$mon]."_".$days[$wday];
		system("mkdir -p ".$work_dir);
		system("mkdir -p ".$tmp_dir);	# recreate tmp directory
	}

	return $work_dir, $tmp_dir;
}

############################################################## 
### to delete

sub read_file {

	my $file = $_[0];
	my $ORFs = {};

	open (F, $file) or die  "Error reading file: $file";
	while (<F>) {
		chomp $_;
		if (/^orf_id/) {next;}
		if (/^ORF_locus/) {next;}
		my @line 			= (split '\t', $_);
#ORF_locus	strand	length	ORF_type	start_codon	reference_annotation	Nterminal_support	dist_from_anno

		my $ORF 		= $line[0];
		my $strand 		= $line[1];
		my $length 		= $line[2];

		my $start_codon = $line[4];
		my $ref 		= $line[5];
		my $nterm		= $line[6];
		my $dist		= $line[7];

		my $ORF_type 	= $line[3];
		#my $ORF_type 	= $line[10];

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
		$ORFs->{$ORF}->{dist} 	= $dist;
	}
	close F; 

	return $ORFs;
}

