#!/usr/bin/perl

use strict;
use warnings;

my %genes;
my %evidence;
while (<>) {
    next if /^\#/o;
    next if /^$/o;
    my @line = split /\t/, $_;
    if ($line[2] eq "gene" && m/ID=([^;]+)/) {
        my $gene = $1;
        $genes{$gene}->{gene} = $_;
    } 
    elsif ($line[2] eq "mRNA" && m/ID=mRNA_([^;]+)/) {
        my $id = $1;
        $genes{$id}->{mRNA} = $_;
        $genes{$id}->{exons} = "";
    }
    elsif ($line[2]=~/CDS|exon/ && m/Parent=mRNA_([^;]+)/){
        my $id = $1;
        $genes{$id}->{exons} = $genes{$id}->{exons} . $_;
        #print "$genes{$id}->{exons}\n\n";
    }
    elsif ($line[2] eq "intron" && $line[1] eq "evidence" && m/Parent=([^;]+)/) {
        #print "$_\n";
        my $id = $1;
        $evidence{$id}++;
    }
}

foreach (keys %evidence) {
    print "$genes{$_}->{gene}" . "$genes{$_}->{mRNA}" . "$genes{$_}->{exons}";
}
