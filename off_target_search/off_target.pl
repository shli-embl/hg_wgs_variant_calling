#!/usr/bin/env perl
use Getopt::Long;
use Bio::SeqIO;
use strict;
my $display_help;
my $input_variants;
my $output_txt;
#my $genome = "/g/steinmetz/shli/data/hg38/hg38.ucsc.fa";
my %genome_seq;
GetOptions('h' => \my $display_help, 'in=s' => \my $input_variants, 'out=s' =>\my $output_txt, 'genome=s' =>\my $genome, 'seq=s' =>\my $grna_seq);
if($display_help)
{
	print "Command: \n\tperl off_target.pl [-h] -in INPUT_VARIANT_FILE -out OUTPUT_FILE -genome GENOME_FASTA -seq GRNA_SEQUENCE\n\n";
	exit 1;
}

if((!$input_variants)||(!$output_txt)||(!$genome)||(!$grna_seq))
{
	print "Missing arguments.\n";
	print "Use -h for looking at usage...\n";
	exit 1;
}


sub mismatch
{
	my $seq1=$_[0];
	my $seq2=$_[1];
	my $length=length($seq1);
	my $no=0;
	for(my $i=0;$i<$length;$i++)
	{
		if(substr($seq1,$i,1) ne substr($seq2,$i,1))
		{
			$no++
		}
	}
	return $no;
}
sub reversecomplement
{
	my $input=$_[0];
	my %rc=("A"=>"T","T"=>"A","C"=>"G","G"=>"C","N"=>"N");

	for(my $p=0;$p<length($input);$p++)
	{
		substr($input,$p,1)=$rc{substr($input,$p,1)};
			
	}
	$input=reverse $input;
	return $input;
}
sub offtarget
{		
	my $seq1=$_[0];
	my $seq2=$_[1];
	my $best_match;
	my $best_match_misno=10000;
	my $best_match_variantpos;
	my $best_match_strand;
	for(my $i = 21;$i<49;$i++)
	{
		if((substr($seq1,$i,2) eq "GG")||(substr($seq1,$i,2) eq "AG"))
		{
			my $off_target_site=substr($seq1,$i-21,20);
			my $tmp_misno=&mismatch($off_target_site,$seq2);
			if($tmp_misno<$best_match_misno)
			{
				$best_match = substr($seq1,$i-21,23);
				$best_match_misno=$tmp_misno;
				$best_match_variantpos=26-$i+21;
				$best_match_strand="+"
			}
		}
	}
	$seq1 = &reversecomplement($seq1);
	for(my $i = 21;$i<49;$i++)
	{
		if((substr($seq1,$i,2) eq "GG")||(substr($seq1,$i,2) eq "AG"))
		{
			my $off_target_site=substr($seq1,$i-21,20);
			my $tmp_misno=&mismatch($off_target_site,$seq2);
			if($tmp_misno<$best_match_misno)
			{
				$best_match = substr($seq1,$i-21,23);
				$best_match_misno=$tmp_misno;
				$best_match_variantpos=26-$i+21;
				$best_match_strand="-"
			}
		}
	}
	return ($best_match,$best_match_misno,$best_match_variantpos,$best_match_strand);
}


print "Start loading genome\n";
my $seqio=Bio::SeqIO->new(-file=>$genome,-format=>"fasta");
while(my $next_seq=$seqio->next_seq)
{
	$genome_seq{$next_seq->display_id}=$next_seq->seq
}
print "Finish loading genome\n";

open(IN,$input_variants);
open(OUT,">".$output_txt);
my $header=<IN>;
chomp $header;
print OUT $header,"\t","subject","\t","best_match","\t","n_of_mismatch","\t","matched_position","\t","matched_strand","\n";
while(my $read_line=<IN>)
{
	chomp $read_line;
	my @read_column=split("\t",$read_line);
	print OUT $read_line,"\t",uc substr($genome_seq{$read_column[0]},$read_column[1]-26,51);
	my @best_match = &offtarget(uc substr($genome_seq{$read_column[0]},$read_column[1]-26,51),$grna_seq);
	print OUT "\t",$best_match[0],"\t",$best_match[1],"\t",$best_match[2],"\t",$best_match[3],"\n";
}
print "Done!\n";
