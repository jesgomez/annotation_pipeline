#!/usr/bin/perl

#use strict;
#use warnings;

my $additional = " ";
my $chunkfile = $ARGV[0];
my $index = $ARGV[1];
my $species = $ARGV[2];
my $alternatives = $ARGV[3];
my $sample = $ARGV[4];
my $gff3 = $ARGV[5];
my $noInFrameStop = $ARGV[6];
my $uniqueGeneId = $ARGV[7];
my $maxtracks = $ARGV[8];
my $strand = $ARGV[9];
my $singlestrand = $ARGV[10];
my $intronlength = $ARGV[11];
my $extrinsic_file = $ARGV[12];
$additional = $ARGV[13];

my $junctions_gff = "junctions.gff";
my $junctions = "junctions.sorted.gff.gz";


my %scaff;
open JUNCTIONS, "<", "$junctions_gff";
while (<JUNCTIONS>){
    chomp;
    my @line = split /\s+/, $_;
    $scaff{$line[0]}++ if (!exists $scaff{$line[0]});
}
close JUNCTIONS;

open FASTA, "<", "$chunkfile";
my $id;
while (<FASTA>) {
    chomp;
    next unless (/^\>/);
    #print "$_\n";
    my @id = split '>', $_;
    $id = $id[1];
    next if (!exists $scaff{$id});
    system "fastafetch -f $chunkfile -i $index -q $id > $id.masked.fa";
    system "tabix $junctions $id > $id.junctions.gff";
    system "augustus --species=$species --hintsfile=$id.junctions.gff --alternatives-from-sampling=$alternatives --sample=$sample --gff3=$gff3 --noInFrameStop=$noInFrameStop --uniqueGeneId=$uniqueGeneId --maxtracks=$maxtracks --strand=$strand --singlestrand=$singlestrand --min_intron_len=$intronlength --extrinsicCfgFile=$extrinsic_file $additional $id.masked.fa > $id.augustus_introns.gff3";
    
}
close FASTA;
