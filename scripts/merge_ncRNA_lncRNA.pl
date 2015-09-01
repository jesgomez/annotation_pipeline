#!/usr/bin/perl

use strict;
use warnings;

my $ncRNA = $ARGV[0];
defined($ncRNA) || die ("##ERROR## This script requires a ncRNA annotation file as first argument\n");

my $lncRNA = $ARGV[1];
defined($lncRNA) || die ("##ERROR## This script requires a lncRNA annotation file as second argument\n");

my $project = $ARGV[2];
defined($project) || die ("##ERROR## This script requires a project name with the assembly version as third argument\n");

my $nc_version = $ARGV[3];
defined($nc_version) || die ("##ERROR## This script requires a version of the nc annotation as last argument\n");


open ncRNA, "<", "$ncRNA";
my %ncRNA;
my %ncRNA_exons;
while (<ncRNA>){
    chomp;
    next if /^\#/o;
    next if /^$/o;
    my @line = split /\t/, $_;
    if ($line[2] eq 'ncRNA') {
        my $id;
        if ($line[8] =~ m/ID=([^;]+)/) {$id = $1;}
        $ncRNA{$line[0]}{$line[3]}{nc}{$id} = $_;
    }
    elsif ($line[2] eq 'exon') {
        my $parent;
        my $id;
        if ($line[8] =~ m/Parent=([^;]+)/) {$parent = $1;}
        if ($line[8] =~ m/ID=/) {$id = $';}
        $ncRNA_exons{$parent}{$id} = $_;
    }
}
close ncRNA;

open lncRNA, "<", "$lncRNA";
my %genes;
my %throw;
my %lnctranscript;
while (<lncRNA>) {
    chomp;
    next if /^\#/o;
    next if /^$/o;
    my @line = split /\t/, $_;   
    if ($line[2] eq 'gene') {
        my $gene;
        if ($line[8] =~ m/ID=/) {$gene = $';}
        $ncRNA{$line[0]}{$line[3]}{lnc}{$gene} = $_;        
    }
    elsif ($line[2] =~ /(ncRNA|transcript)/i) {
        my $id;
        if ($line[8] =~ m/ID=([^;]+)/) {$id = $1;}
        my @gene = split /T/, $id;
        if ($line[8] =~ m/small_ncRNA/) {$throw{$gene[0]}++;}
        $lnctranscript{$gene[0]}{$id} = $_;
    }
    elsif ($line[2] eq 'exon') {
        my $parent;
        my $id;
        if ($line[8] =~ m/Parent=([^;]+)/) {$parent = $1;}
        if ($line[8] =~ m/ID=([^;]+)/) {$id = $1;}
        $ncRNA_exons{$parent}{$id} = $_;
    }

}
close lncRNA;

print "##gff-version 3\n";
my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst)=localtime(time);
printf "# date: %4d-%02d-%02d %02d:%02d:%02d\n",$year+1900,$mon+1,$mday,$hour,$min,$sec;

my $i = 1;
foreach my $scaffold (sort keys %ncRNA) {
    my $gene_id;
    foreach my $start (sort {$a<=>$b} keys %{$ncRNA{$scaffold}}) {
        if (exists $ncRNA{$scaffold}{$start}{nc}) { 
            foreach my $transcript (keys %{$ncRNA{$scaffold}{$start}{nc}}) {
                $gene_id = $project . "nc" . $nc_version . sprintf("%06d",$i);
                my @line = split /\t/, $ncRNA{$scaffold}{$start}{nc}{$transcript};
                my $type;
                if ($line[8] =~ m/Type=/) {$type = $';}
                print "\#\#\#\n$scaffold\tCNAG\tgene\t$start\t$line[4]\t.\t$line[6]\t$line[7]\tID=$gene_id\n";
                print "$scaffold\tCNAG\ttranscript\t$start\t$line[4]\t$line[5]\t$line[6]\t$line[7]\tID=$gene_id" . "T1;Name=$gene_id" . "T1;Description=$type\n";
                foreach my $exons (sort keys %{$ncRNA_exons{$transcript}}) {
                    my @line = split /\t/, $ncRNA_exons{$transcript}{$exons};
                    my $exon;
                    if ($line[8] =~ m/.exon/) {$exon = $';}
                    print "$scaffold\tCNAG\texon\t$line[3]\t$line[4]\t$line[5]\t$line[6]\t$line[7]\tParent=$gene_id" . "T1;ID=$gene_id" . "T1.exon" . "$exon;Name=$gene_id" . "T1\n";
                }
            }
            $i++;
        }
        if (exists $ncRNA{$scaffold}{$start}{lnc}) {
            foreach my $gene (keys %{$ncRNA{$scaffold}{$start}{lnc}}) {
                if (!exists $throw{$gene}) {
                    $gene_id = $project . "nc" . $nc_version . sprintf("%06d",$i);  
                    my @line = split /\t/, $ncRNA{$scaffold}{$start}{lnc}{$gene};
                    $line[8] = "ID=$gene_id";
                    local $" = "\t";
                    print "\#\#\#\n@line\n";  
                    foreach my $transcript (sort keys %{$lnctranscript{$gene}}) { 
                        my $t;
                        my $old_id;
                        my $attributes;
                        my @line = split /\t/, $lnctranscript{$gene}{$transcript};
                        if ($line[8] =~ m/Name=([^;]+)/) {$old_id = $1; $attributes = $'}                 
                        if ($old_id =~ m/T/) {$t = $';}              
                        print "$scaffold\tCNAG\ttranscript\t$line[3]\t$line[4]\t$line[5]\t$line[6]\t$line[7]\tID=$gene_id" . "T" . "$t;Name=$gene_id" . "T" . "$t" . "$attributes\n";
                        foreach my $exons (sort keys %{$ncRNA_exons{$transcript}}) {
                            my @line = split /\t/, $ncRNA_exons{$transcript}{$exons};
                            my $exon;
                            if ($line[8] =~ m/.exon([^;]+)/) {$exon = $1;}
                            print "$scaffold\tCNAG\texon\t$line[3]\t$line[4]\t$line[5]\t$line[6]\t$line[7]\tParent=$gene_id" . "T" . "$t;ID=$gene_id" . "T" . "$t.exon" . "$exon;Name=$gene_id" . "T" . "$t\n";
                        }
                    }
                }
                $i++;
            }
        }
    }
}
