#!/usr/bin/perl

use strict;
use warnings;

my $fastafile = $ARGV[0];
#defined($fastafile) || die ("##ERROR## This script requires a fasta file as argument\n");

my $dir = $ARGV[1]; # A directory where fasta files will be generated

my $trained = $ARGV[2];


open INFILE, "<", "$fastafile"; #|| die ("##ERROR## Cannot open file $fastafile\n");

my $i = -1;

while (<INFILE>) {
   chomp;
   if ($_ =~ />/) {
      close OUT;
      $i++;
      my @scaffold = split //, $_;
      shift @scaffold;
      my $scaffold = join '', @scaffold;
      open OUT, ">", "$dir/scaffold_$i.fa";

      print OUT ">$scaffold\n";
   }
   else {
      print OUT "$_\n";
  
   }
}

close OUT;

for (my $n=0; $n <= $i; $n++) {
   system "/project/devel/aateam/src/GlimmerHMM/bin/glimmerhmm_linux_x86_64 '$dir'/scaffold_'$n'.fa -d $trained -g > '$dir'/scaffold_'$n'.gff3;"
}
