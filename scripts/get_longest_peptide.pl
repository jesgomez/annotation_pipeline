#!/usr/bin/env perl

############################################################
# script to get the isoform of a gene that codes for the longest peptide.
#
############################################################

use strict;
my $file = $ARGV[0];
#now lets get the longest protein per gene
open PIN,"cat $file | FastaToTbl | ";
open OUT,"|TblToFasta";
my %seq;
my %pid;
while (<PIN>){
  my ($id,$seq)=split;
  $id=~/(\S+)[TP]\d+/;
  if (exists $seq{$1}){
    if (length($seq)>length($seq{$1})){
      $pid{$1}=$id; 
    }
  }else{
    $seq{$1}=$seq;
    $pid{$1}=$id; 
  }
}
foreach my $g (sort keys %seq){
  print OUT $pid{$g},"\t",$seq{$g},"\n";
}
close PIN;
close OUT;

