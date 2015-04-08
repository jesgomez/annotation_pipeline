#!/usr/bin/perl

use strict;
use warnings;

my $chunkfile = $ARGV[0];
my $index = $ARGV[1];
my $junctions_gff = $ARGV[2];
my $junctions = $ARGV[3];
my $geneid_parameters = $ARGV[4];
my $geneid_options = $ARGV[5];

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
    #print "$id\n";
    system "fastafetch -f $chunkfile -i $index -q $id > $id.masked.fa";
    system "tabix $junctions $id > $id.junctions.gff";
    system "geneid -R $id.junctions.gff -P $geneid_parameters $geneid_options $id.masked.fa > $id.geneid_with_introns.gff3"; 
}
close FASTA;
