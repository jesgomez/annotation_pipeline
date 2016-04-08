#! /usr/bin/perl

use warnings;
use strict;

my $gff = $ARGV[0];
open GFF, "<", "$gff";
my %gene;
my $id;
my $start;
my $scaffold; 
print STDERR "parsing gff\n";
while (<GFF>){
    chomp;
   # print "$_\n";
    my @line = split /\t/, $_;
    if ($line[2]=~ m/gene/){
        if (~m/locus_tag \"([^\"]+)/){
            $id = $1;
        }
        $start = $line[3];
        $scaffold=$line[0];
        $gene{$scaffold}{$start}{$id}="$_";
        print STDERR "$id\n";
        #print "$_\n";
    }
    elsif ($line[2]=~ m/repeat_region/){
        if (~m/locus_tag \"([^\"]+)/){
            $id = $1;
        }
        my $line = join "\t", @line;
        $start = $line[3];
        $scaffold=$line[0];
        $gene{$scaffold}{$start}{$id}="$line";
        print STDERR "$id\n";
        #print "$_\n";
    }
    else {
       # print "$_\n" if (!exists $gene{$scaffold}{$start}{$id});
        $gene{$scaffold}{$start}{$id}="$gene{$scaffold}{$start}{$id}" . "\n" . "$_";
    }
  
}
print STDERR "sorting\n";

foreach my $scaff (sort keys %gene){
     foreach my $coord (sort {$a<=>$b} keys %{$gene{$scaff}}){
         foreach my $loc (keys %{$gene{$scaff}{$coord}}){
             #print "$scaff\t$coord\t$loc\n";
             print "$gene{$scaff}{$coord}{$loc}\n";
         }
     }
}
