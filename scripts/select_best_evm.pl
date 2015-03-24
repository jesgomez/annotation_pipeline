#!/usr/bin/perl

use strict;

my $gmap_alignments = $ARGV[0];
#print "$gmap_alignments\n";

my $total_evm = `ls weights_*.txt | wc -l`;
chomp $total_evm;
my %bp;
my $i;
my @line;


for ($i=1; $i<=$total_evm; $i++) {
   #print "$i\n";
   system "intersectBed -wo -a evm_weights_'$i'.TEcleaned.CDS.gff3 -b ../'$gmap_alignments'  > evm_weights_'$i'.TEcleaned.BT.gmap.out";
   open BEDTOOLS, "<", "evm_weights_$i.TEcleaned.BT.gmap.out" ;
   $bp{$i} = 0;
   while (<BEDTOOLS>) {
      chomp;
      @line = split /\t/, $_;
      my $current = pop @line;
     # print "$current\n";
      $bp{$i} = $bp{$i} + $current;
    #  print "$bp{$i}\n";
   }
  # print "$bp{$i}\n";
   close BEDTOOLS; 
 
}

my $best = 0;
my $bestid = "";
foreach (keys %bp) {
  # print "$_\n";
  # print "$bp{$_}\n";
  # print "$best\n";
   if ($bp{$_} > $best) {
      $best = $bp{$_};
      $bestid = $_;
   }
}

system "ln -s evm_weights_'$bestid'.TEcleaned.gff3 evm_weights_'$bestid'.best.TEcleaned.gff3";
print "evm_weights_$bestid.TEcleaned.gff3\n";
