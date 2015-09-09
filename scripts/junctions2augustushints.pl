#!/usr/bin/perl

############################################################
#
# script to make the junctions gff3 to the format accepted by Augustus.
#
# Author: J. Gómez Garrido
############################################################

use strict;
use warnings;

#my $junctionfile = shift @ARGV;
#defined($junctionfile) || die ("##ERROR## This script requires a gff3 file as argument\n");

#open JUNCTIONS, "<", "$junctionfile" || die ("##ERROR## Cannot open file $junctionfile\n");

while (<STDIN>) {
    chomp;
    next if /^\#/o;
    my @line = split /\t/, $_;
   # my @mult = split /\./, $line[5];
   # my $mult = $mult[0];
   # print "$mult\n";
    print "$line[0]\t$line[1]\t$line[2]\t$line[3]\t$line[4]\t$line[5]\t$line[6]\t$line[7]\tsrc=E\n";
}

