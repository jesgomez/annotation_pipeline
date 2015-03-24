#!/usr/bin/env perl
use lib "~talioto/myperlmods/bioperl/bioperl-live/";
use Bio::Tools::GFF;
use Getopt::Long;
#use GFFfeat;
use Bio::SeqIO;
use Bio::Seq;
use Bio::FeatureIO;
use Bio::FeatureIO::gff;
use Bio::SeqFeature::Annotated;
use Bio::SeqFeature::Gene::Transcript;
use Bio::SeqFeature::Gene::Exon;
use strict;
#can not run on exonpairs... must run on full gene predictions.
# MUST be GFF2 format!
# gff data structure
# my $seqname  = 0;
# my $source  = 1;
# my $feature  = 2;
# my $start    = 3;
# my $end      = 4;
# my $score    = 5;
# my $strand   = 6;
# my $frame    = 7;
# my $group    = 8;
my $SCORE = 0;
my $GROUP = '';
my $DONOR = '';
my $FRAME = ".";
my $donor = 0;
my $STRAND = "+";
my $SEQNAME = "";
#my @gff = ();
my $u12isdonor = 0;
my $feature_type = 'exon';
my $version = 3;
my $oversion = 3;
my $userscore = undef;
my $minintronlength = 15;
GetOptions(
	   'v:s'       => \$version,
	   'ov:s'      =>$oversion,
	   'us:s'      => \$userscore,
	   'f:s'       => \$feature_type
	  );
#
#$oversion = $version unless $oversion;
print STDERR "Reading GFF input...\n";
my $gffio = Bio::Tools::GFF->new(-fh => \*STDIN, -gff_version => $version); 
#$gffio->close();
my $gffout = Bio::Tools::GFF->new(-fh => \*STDOUT, -gff_version => $oversion); 
my $gfferror = Bio::Tools::GFF->new(-fh => \*STDERR, -gff_version => $oversion); # loop over the input stream
#my %transcripts;
my $current_transcript = Bio::SeqFeature::Gene::Transcript->new();
my $current_transcript_name = undef;
while (my $feature = $gffio->next_feature()) {
    next if ($feature->primary_tag !~ /$feature_type/); # CDS|exon|First|Internal|Terminal|Single|HSP|match_part
    #print STDERR $feat->seq_id, "\n";
    #next if $feature->seq_id eq "ChrUn.10757";
#     foreach my $tag ($feature->get_all_tags()){
#     print STDERR $tag,"\n";
# }
    my %att;
    my @fields = split /\t/,$feature->gff_string;
    foreach my $key ( $feature->get_all_tags ) {
	my @values = $feature->get_tag_values($key);
	$att{$key}=\@values;
    }
 #    my $ac = $feature->annotation();
#     my %att;
    #my $stcf = 0;
    #print STDERR $feature->gff_string,"\n";
#     my @fields = split /\t/,$feature->gff_string;
#     foreach my $key ( $ac->get_all_annotation_keys() ) {
# 	#print STDERR "KEY:$key\n";
# 	my @values = $ac->get_Annotations($key);
# 	$att{$key}=\@values;
#     }
    #      if ($feature->primary_tag eq "mRNA"){
    #        $transcripts{$att{ID}}->{mRNA}=$feature;
    #      }
    if (exists $att{'Parent'}) {
	$att{transcript_id}=$att{Parent};
    } elsif (exists $att{transcript_id}) {
	#print STDERR "TID:".$att{transcript_id}->[0]."\n";
    } elsif (exists $att{group}) {
	$att{transcript_id}=$att{group};
	#print STDERR "GROUP:".$att{transcript_id}->[0]."\n";
    } elsif (exists $att{GenePrediction}) {
	$att{transcript_id}=$att{GenePrediction};
	#print STDERR "GROUP:".$att{transcript_id}->[0]."\n";
    } elsif (exists $att{Match}) {
	$att{transcript_id}=$att{Match};
	#print STDERR "GROUP:".$att{transcript_id}->[0]."\n";
    } else {
	foreach my $k (keys %att){print STDERR "$k=",join(",",@{$att{$k}}),"\n";}
	my $grp = $fields[8];
	$grp =~ s/\s+/_/g;
	$att{transcript_id}=$grp;
    }
    foreach my $p (@{$att{transcript_id}}) {
	#print STDERR "Parent:$p\n";
	my $newfeature = Bio::SeqFeature::Gene::Exon->new(-start=>$feature->start,
							  -end=>$feature->end,
							  -seq_id=>$feature->seq_id,
							  -source_tag=>$feature->source_tag,
							  -primary_tag=>$feature->primary_tag,
							  -strand=>$feature->strand,
							  -frame=>$feature->frame
							 );
	# my $original_p = $p;
	# $p = $feature->seq_id ."_". $feature->strand() ."_$p";
	foreach my $key ( $feature->get_all_tags ) {
	    my @values = $feature->get_tag_values($key);
	    $att{$key}=\@values;
	    foreach my $val (@values) {
		$newfeature->add_tag_value($key,$val);
		# if ($key eq "Parent") {
		#     $newfeature->add_tag_value($key,$feature->seq_id ."_". $feature->strand() ."_$val");
		# } else {
		#     $newfeature->add_tag_value($key,$val);
		# }
	    }
	}
# 	$p = $feature->seq_id ."_". $feature->strand() ."_$p";
# 	foreach my $key ( $ac->get_all_annotation_keys() ) {
# 	    my @values = $ac->get_Annotations($key);
# 	    foreach my $val (@values) {
# 		#print STDERR "$val\n";
# 		if ($key eq "Parent") {
# 		    $newfeature->add_tag_value($key,$feature->seq_id ."_". $feature->strand() ."_$val");
# 		} else {
# 		    $newfeature->add_tag_value($key,$val);
# 		}
# 	    }
# 	}

	if ((defined $current_transcript_name) && ($current_transcript_name eq $p)) {
	    $current_transcript->add_exon($newfeature);
	} else {
	    printLastTranscript($current_transcript) if defined $current_transcript_name;
	    $current_transcript = undef;
	    $current_transcript = Bio::SeqFeature::Gene::Transcript->new();
	    $current_transcript->add_exon($newfeature);
	    $current_transcript_name = $p;
	}
    }
}
$gffio->close;
printLastTranscript($current_transcript);

sub printLastTranscript{
    my $t = shift;
    my @lastGRPgff;
    #print STDERR "Transcript:$t\n";
    my $exon_number = 0;
    my $tid = "";
    my $tlength = 0;
    foreach my $e (sort {$a->start <=> $b->start} $t->exons) {$tlength+=$e->length;}
    #print STDERR "Transcript is $tlength bp long\n";
    my $partial_length = 0;
    foreach my $feature (sort {$a->start <=> $b->start} $t->exons) {
	$exon_number++;
	#next if $feature->seq_id eq "ChrUn.10757";
	#print STDERR $feature->gff_string,"\n";
	my $ac = $feature->annotation();
	my %att;
	foreach my $key ( $feature->get_all_tags) {
	    my @values = $feature->get_tag_values($key);
	    $att{$key}=\@values;
	}
# 	foreach my $key ( $ac->get_all_annotation_keys() ) {
# 	    #print STDERR "KEY:$key\n";
# 	    my @values = $ac->get_Annotations($key);
# 	    $att{$key}=\@values;
# 	}
	#      if ($feature->primary_tag eq "mRNA"){
	#        $transcripts{$att{ID}}->{mRNA}=$feature;
	#      }
	if (exists $att{'Parent'}) {
	    $att{transcript_id}=$att{Parent};
	} elsif (exists $att{transcript_id}) {
	    #print STDERR "TID:".$att{transcript_id}->[0]."\n";
	} elsif (exists $att{group}) {
	    $att{transcript_id}=$att{group};
	    #print STDERR "GROUP:".$att{transcript_id}->[0]."\n";
	} elsif (exists $att{GenePrediction}) {
	    $att{transcript_id}=$att{GenePrediction};
	    #print STDERR "GROUP:".$att{transcript_id}->[0]."\n";
	} elsif (exists $att{Match}) {
	    $att{transcript_id}=$att{Match};
	    #print STDERR "GROUP:".$att{transcript_id}->[0]."\n";
	} else {
	    #foreach my $k (keys %att){print STDERR "$k=",join(",",@{$att{$k}}),"\n";}exit;
	    #my $grp = $fields[8];
	    #$grp =~ s/\s+/_/g;
	    #$att{transcript_id}=$grp;
	}

	$tid = $att{transcript_id}->[0];
	if ($exon_number == 1) {
	    $GROUP = $att{transcript_id}->[0];
	    $STRAND = $feature->strand;
	    $SEQNAME = $feature->seq_id;

	    $DONOR=$feature->end +1;
	    $SCORE = $feature->score if defined $feature->score;
	    #print STDERR "exon 1 is ".$feature->length." bp long\n";
	    $partial_length += $feature->length;

	} else {
	    $donor = $feature->end +1;
	    my $newintron = Bio::SeqFeature::Generic->new(-start=>$DONOR,
							  -end=>$feature->start -1,
							  -seq_id=>$feature->seq_id,
							  -source_tag=>$feature->source_tag,
							  -primary_tag=>"Intron",
							  -strand=>$feature->strand,
							  -score=>$feature->score,
							  -frame=>$feature->frame,
							  #-annotation=>$e->annotation
							 );
	    my $has_target = 0;
	    foreach my $annot (keys %att) {
		foreach my $val (@{$att{$annot}}) {
		    #print STDERR "VALUE:$val\n";
		    if ($annot eq "Target"){
			$has_target=1;
			$annot="target";#$val =~s/\s+$//;
			
			my ($n,$st,$en,$tstrand)= split /\s+/,$val;
			#print STDERR "$val\n".$feature->strand."\n$tstrand\n";
			if ((($feature->strand > 0)&&($tstrand eq "+"))||(($feature->strand < 0)&&($tstrand eq "-"))){
			    #$st = $en; $en = $st+1;
			    $en = $st; $st = $st-1;
			}else{
			    #$en = $st; $st = $st-1;
			    $st = $en; $en = $st+1;
			}
			$val = "$n $st $en $tstrand";
		    }
		    $newintron->add_tag_value($annot,$val);
		}
	    }
	    if (!$has_target){
		if ($feature->strand > 0){
		    $newintron->add_tag_value('target',$tid." ".$partial_length." ".($partial_length + 1)." +");
		}else{
		    $newintron->add_tag_value('target',$tid." ".($tlength - $partial_length)." ". ($tlength - $partial_length+1)." +");
		}
	    }


	    $newintron->end($feature->start -1);
	    $newintron->start($DONOR);
	    $DONOR=$donor;
	    $partial_length += $feature->length;
	    #print STDERR "exon is ".$feature->length." bp long\n";
	    $newintron->score(($SCORE + $feature->score)/2) if defined $feature->score;
	    $SCORE = $newintron->score if defined $newintron->score;
	    if ($feature->frame ne ".") {
		if ($newintron->strand == 1) {
		    $newintron->frame((3 - $feature->frame)%3 );
		} else {
		    $newintron->frame((($feature->end - $feature->start + 1) - $feature->frame)%3 );
		}
	    }
	    $newintron->score($userscore) if ! defined $newintron->score  && defined $userscore;
	    #print STDERR $newintron->gff_string,"\n";
	    if (( ($newintron->end - $newintron->start + 1) > 5000000 )||( ($newintron->end - $newintron->start + 1) < $minintronlength )) {
		print STDERR ('Warning: abnormal intron size (',($newintron->end - $newintron->start + 1),")\n");
		$gfferror->write_feature($newintron);
	    }else{
		push @lastGRPgff, $newintron;
	    }
	}

    }
    if (scalar @lastGRPgff) {

	if ($lastGRPgff[0]->strand eq "1") {
	    my $counter = 1;
	    foreach my $int (@lastGRPgff) {
		$int->remove_tag('ID');
		$int->add_tag_value('ID',$tid."_intron_$counter");
		$int->add_tag_value('inum',$counter);
		$gffout->write_feature($int);
		$counter++;
	    }
	} else {
	    my $counter = scalar @lastGRPgff;
	    foreach my $int (@lastGRPgff) {
		$int->remove_tag('ID');
		$int->add_tag_value('ID',$tid."_intron_$counter");
		$int->add_tag_value('inum',$counter);
		$gffout->write_feature($int);
		$counter--;
	    }
	}
	@lastGRPgff = ();
    }








}

#$gffio->close();
#   if (scalar @lastGRPgff) {

#     if ($lastGRPgff[0]->strand eq "1") {
#       my $counter = 1;
#       foreach my $int (@lastGRPgff) {
# 	$int->add_tag_value('inum',$counter);
# 	$gffout->write_feature($int);
# 	$counter++;
#       }
#     } else {
#       my $counter = scalar @lastGRPgff;
#       foreach my $int (@lastGRPgff) {
# 	$int->add_tag_value('inum',$counter);
# 	$gffout->write_feature($int);
# 	$counter--;
#       }
#     }
#   }
$gffout->close;
$gfferror->close;



