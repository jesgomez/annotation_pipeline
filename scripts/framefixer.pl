#!/usr/bin/env perl
use strict;
use Getopt::Long;
use lib "/project/devel/aateam/perlmods";
use lib "/apps/BIOPERL/1.6.1/lib/perl5";
use Bio::Seq;
use Bio::SeqIO;
use Bio::DB::Fasta;
use File::Basename qw( fileparse );
use SeqOp;

# Need to add GFF2 FITS (geneid-like) output
my $version = 3;
my $source = 0;
my $removestop = 0;
my $addstop = 0;
my $txtag = 'mRNA';
my $mrna = 0;
my $ov = 3;
my $fits = 0;
my $fl = 0;
my $addutr = 0;
my $addname = 1;
my $printgene = 1;
my $nr = 0;
my $fa= 0;
my $code = 1;
my $verbose = 0;
my $geneid_param = "/home/devel/talioto/param/human.090903.param";
GetOptions(
	   'v|version:s'      => \$version,
	   'rmstop|nostop'  => \$removestop,
	   'addstop!'        => \$addstop,
	   'UTR'        => \$addutr,
	   'addname!'        => \$addname,
	   'gene!'           => \$printgene,
	   'txtag:s'        => \$txtag,
	   'mrna'           => \$mrna,
	   'ov:s'           => \$ov,
	   'fits'           => \$fits,
	   'fl'             => \$fl,
	   'nr'             => \$nr,
           's|seq:s'          => \$fa,
           'c|code:s'           => \$code,
	   'v'              => \$verbose,
	   'param:s'        => \$geneid_param
	   );

# open the first file: table of codon usage frequencies

open(PARAM,"<$geneid_param") or die "$0: the file $geneid_param can not be opened: $!\n";

# load the frequencies of codon usage into the hash table %pcodons
# this hash is indexed by a triplet of nucleotides or codon

my %llhex;
while (my $line = <PARAM>) {
  next if $line!~/Markov_Transition_probability_matrix/;
  while(my $entry = <PARAM>){
    last if $entry!~/\w/;
    last if $entry=~/^#/;
    chomp $entry;
    my($hex,$i,$frame,$lls)=split /\s+/,$entry;
    $llhex{$hex}{$frame}=$lls;
  }     
}
close(PARAM);

my $db = Bio::DB::Fasta->new($fa);
if ($mrna){$txtag = 'mRNA';}
my %transcripts;
my %transcript_exons;
my %tx_gene;
my %genes;
my %tx;
my %tx_target;
my %printed;
my $check_start_stop_codon = 0;
#if($addutr){print STDERR "will at UTRs if present\n";}
while(<>){
    next if m/^[# ]/;
    chomp;
    my @gtf = split "\t",$_;
    next if $gtf[2] !~/CDS|First|Internal|Terminal|initial|exon|gene|transcript|RNA|start_codon|stop_codon/i;
    my %att;
    if ($version eq "2.5"){
	#print STDERR $gtf[8],"\n";
	$gtf[8]=~s/;\s*$//;
	my @list = split ";",$gtf[8];
	foreach my $kv (@list) {
	    $kv =~ /(\S+)\s+(\S.*);?/;
	    $att{$1}=$2;
	}
	$gtf[8]='';
	
	if (exists $att{transcript_id} && exists $att{Parent}){delete $att{Parent};}
	foreach my $a (sort keys %att) {
	    $att{$a}=~s/;$//;
	    $att{$a}=~s/^\s*//;
	    $att{$a}=~s/"$//;
	    $att{$a}=~s/^"//;
	    #$gtf[8].="$a=$att{$a};";
	    if ($gtf[2]=~/gene/){
	      if ($a eq 'gene_id') {
		$gtf[8].="ID=$att{$a};";
	      }else{
		$gtf[8].="$a=$att{$a};";
	      }
	      
	    }elsif ($gtf[2]=~/transcript|mRNA/){
	      if ($a eq 'transcript_id') {
		$gtf[8].="ID=$att{$a};";
	      }elsif ($a eq 'gene_id') {
		$gtf[8].="Parent=$att{$a};";
	      }else{
		$gtf[8].="$a=$att{$a};";
	      }
	      
	    }else{
	      if ($a eq 'transcript_id') {
		$gtf[8].="Parent=$att{$a};";
	      }elsif ($a eq 'gene_id') {
		#$gtf[8].="gene=$att{$a};";
	      }else{
		$gtf[8].="$a=$att{$a};";
	      }
	  }
	}
	$gtf[8]=~s/;$//;
	if ($gtf[2] eq "gene" && (exists $att{gene_id} )){
	  $gtf[8]="ID=$att{gene_id}";
	    $genes{$att{gene_id}}=\@gtf;
	}elsif($gtf[2] =~/transcript|mRNA/ && (exists $att{gene_id} )){
	    $tx_gene{$att{transcript_id}}=$att{gene_id};
	}elsif($gtf[2] =~/codon/ && (exists $att{gene_id} )){
	    $check_start_stop_codon = 1;
	}

    }elsif ($version eq "2"){
	#print STDERR $gtf[8],"\n";
	my $grp = $gtf[8];
	$gtf[8]="Parent=$gtf[8];";
	$att{transcript_id}=$grp;
	if ($gtf[2] eq "gene"){
	    $gtf[8]="ID=$grp;";
	    $genes{$grp}=\@gtf;
	}elsif ($gtf[2] =~/First|Internal|Terminal|Single|initial|internal|terminal/i){
	    $gtf[2] = 'CDS';
	}
	
    }elsif($version eq "3"){
	#print STDERR $gtf[8],"\n";
	if ($gtf[2] eq "gene"){
	    if ($gtf[8]=~m/ID=([^;]+)/){
		$genes{$1}=\@gtf;
	    }
	}elsif($gtf[8]=~m/Parent=([^;]+)/){
	    $att{transcript_id}=$1;
	    #print "$1\n";
	}
    }
    if ($gtf[2] =~/transcript|RNA/){
	$gtf[2] =$txtag;
	my $gene =0;
	my $id =0;
	if ($gtf[8]=~m/Parent=([^;]+)/){
	    $gene = $1;
	}
	if ($gtf[8]=~m/ID=([^;]+)/){
	    $id = $1;
	}
	if ($id && $gene){
	    $tx_gene{$id}=$gene;
	    #print STDERR $gtf[7],"\n";
	    $tx{$id}=\@gtf;
	    #print STDERR "$id\t$gene\n";
	}
    }elsif ($gtf[2] =~/CDS|First|Internal|Terminal|initial|exon|start_codon|stop_codon/i && exists $att{transcript_id}) {
      if ($gtf[2] =~/CDS|exon/){
	$gtf[8] =~s/ID=([^;]+);?//;
	push @{$transcript_exons{$att{transcript_id}}->{exons}},\@gtf;
      }
      if ($gtf[8] =~/Target=([^;\s]+)/){
	 $tx_target{$att{transcript_id}}=$1;
      }

      push @{$transcripts{$att{transcript_id}}},\@gtf;
    }

}
my @unsorted;


if ($nr) {
  print STDERR "\nRemoving redundant transcripts...";
  my %unique_transcripts;
  foreach my $t (sort keys %transcript_exons) {
    my @exons = sort {$a->[3] <=> $b->[3]} @{$transcript_exons{$t}->{exons}};
    my $location_string = "";
    foreach my $exon (@exons) {
      $location_string .= ($exon->[0] ."_". $exon->[3] ."_". $exon->[4] ."_". $exon->[6] ."_");
    }
    # print STDERR  "$t\t$location_string\n";
    if (! exists $unique_transcripts{$location_string}) {
      $unique_transcripts{$location_string}=$t;
    }
  }
  foreach my $ls (keys %unique_transcripts) {
    #my @CDS = sort {$a->start <=> $b->start} @{$t->{CDS}};
    #$t->{CDS}=\@CDS;
    push @unsorted, $transcripts{$unique_transcripts{$ls}};
    print STDERR "$unique_transcripts{$ls}\n";
  }
} else {
  foreach my $t (keys %transcripts) {
    #my @CDS = sort {$a->start <=> $b->start} @{$t->{CDS}};
    #$t->{CDS}=\@CDS;
    push @unsorted, $transcripts{$t};
  }
}


my @sorted_transcripts = sort gffsort @unsorted;


foreach my $transcript (@sorted_transcripts){
    my @tgff = @{$transcript->[0]};
    $tgff[2]= $txtag;
    $tgff[4]= $transcript->[-1]->[4];
    $tgff[5]='.';
    $tgff[8]=~/ID=([^;]*)/;
    my $originalid = $1;
    #print STDERR "original $originalid\n";
    $tgff[8]=~/Parent=([^;]*)/;
    my $tid = $1;
    $tgff[8]="ID=$1";

    $tgff[7]='.';
    if ($ov == 2){$tgff[8]=$1;}
   
    if ($ov !=2){
      if (exists $tx_gene{$tid}) {
	if (exists $genes{$tx_gene{$tid}}) {
	  if ($printgene){
	    if(! (exists ($printed{$tx_gene{$tid}}))){
	      print (join("\t",@{$genes{$tx_gene{$tid}}}),"\n");
	      $tgff[8].=";Parent=$tx_gene{$tid}";
	      $printed{$tx_gene{$tid}}++;
	    }
	  }
	  
	} elsif(exists $tx_gene{$tid}){
	  if(! (exists ($printed{$tx_gene{$tid}}))){
	    print (join("\t",($tgff[0],$tgff[1],'gene',$tgff[3],$tgff[4],'.',$tgff[6],'.',"ID=".$tx_gene{$tid})),"\n");
	    $printed{$tx_gene{$tid}}++;
	  }
	}else {
	  if ($printgene){	  
	    $tgff[8].=";Parent=g_$tid";
	    print (join("\t",($tgff[0],$tgff[1],'gene',$tgff[3],$tgff[4],'.',$tgff[6],'.',"ID=g_$tid")),"\n");
	  }
	}
	delete $tx_gene{$tid};
      } else {
	    #print STDERR "HERE\n";
	if ($printgene){
	  $tgff[8].=";Parent=g_$tid";
	  print (join("\t",($tgff[0],$tgff[1],'gene',$tgff[3],$tgff[4],'.',$tgff[6],'.',"ID=g_$tid")),"\n");
	}

      }
    }
    my @cds;
    my @exons;
    my $count = 1;
    my $numexon = 0;
    my $numcds = 0;
    
    ### Make sure CDS starts and ends with stop. If not it's a partial annotation/prediction and will be call Internal
    my $stop = 0;
    my $start = 0;
    my $numstops = 0;
    my $numstarts = 0;
    foreach my $gffrec (@{$transcript}){
	$numexon++ if $gffrec->[2] eq "exon";
	$numcds++ if $gffrec->[2] eq "CDS";
	if ($gffrec->[2] eq "start_codon"){
	  #print STDERR "start $gffrec->[2]\n";
	  $numstarts++; 
	  if ($gffrec->[6] eq '+'){
	    $start = $gffrec->[3];
	  }else{
	    $start = $gffrec->[4];
	  }
	}
	if ($gffrec->[2] eq "stop_codon"){
	  $numstops++; 
	  if ($gffrec->[6] eq '+'){
	    $stop = $gffrec->[4];
	  }else{
	    $stop = $gffrec->[3];
	  }
	}
    }
    
    foreach my $gffrec (@{$transcript}){
      $gffrec->[8]=~s/;$//;
	if ($gffrec->[3] < $tgff[3]){
	  $tgff[3]=$gffrec->[3];
				     #print STDERR "Changed tx coord\n";
	}
	if ($gffrec->[4] > $tgff[4]){
	  $tgff[4]=$gffrec->[4];#print STDERR "Changed tx coord\n";
	}
	if ($gffrec->[2] eq "CDS"){
	  if($ov == 2){
	    $gffrec->[8]="$tid";
	  }else{
	    $gffrec->[8].=";ID=$tid"."_cds";
	  }
	    push @cds, $gffrec;
	}elsif($gffrec->[2] eq "exon"){
	  if($ov == 2){
	    $gffrec->[8]="$tid";
	  }else{
	    if ($gffrec->[6] eq "+"){
		$gffrec->[8].=";ID=$tid"."_exon$count";
	    }else{
		$gffrec->[8].=";ID=$tid"."_exon".($numexon-$count + 1);
	    }
	    $count++;
	  }
	    push @exons, $gffrec;
	}else{
	  $gffrec->[7]='.';
	}
    }
    if($ov != 2){
      if (exists $tx{$tid}){
	if($addname){
	  my $name = $tid;
	  if(exists $tx_target{$tid}){
	    $name = $tx_target{$tid};
	  }
	  if($tx{$tid}->[8] =~/Name=/){
        	$tx{$tid}->[8] =~s/Name=[^;]+/Name=$name/;
	  }else{
	  	$tx{$tid}->[8].=";Name=$name";
	  }
	}
	print (join("\t",@{$tx{$tid}}),"\n");
      }else{
	if($addname){
	  my $name = $tid;
	  if(exists $tx_target{$tid}){
	    $name = $tx_target{$tid};
	  }
	  $tgff[8].=";Name=$name";
	}
	print (join("\t",@tgff),"\n");
      }
    }

    my $cdsnew= cdsPhase($db, \@cds);
    my @scds = sort gffsort @$cdsnew; #make sure cds records are sorted
    my @sexons = sort gffsort @exons; #make sure exon records are sorted
    my $CDSstart = $scds[0]->[3];
    my $CDSend =  $scds[-1]->[4];
    
    ### Check for stop

    if ($check_start_stop_codon) {
      next if $fl && !($numstarts && $numstops);
      if ($addstop && @scds && $numstops == 1) {
	if ($tgff[6] eq "+") {
	  $scds[-1]->[4]+=3;
	} else {
	  $scds[0]->[3]-=3;
	}
      }
      if ($fits) {
	if ($numcds == 1) {
	  if ($numstarts && $numstops) {
	    $scds[0]->[2] = 'Single';
	  } elsif ($numstarts) {
	    $scds[0]->[2] = 'First';
	  } elsif ($numstops) {
	    $scds[0]->[2] = 'Terminal';
	  } else {
	    $scds[0]->[2] = 'Internal';
	  }
	} else {
	  for (my $i = 0;$i<$numcds;$i++) {
	    $scds[$i]->[2] = 'Internal'; 
	    if (($i == 0)&&($scds[$i]->[6] eq '+')) {
	      $scds[$i]->[2] = 'First' if ($numstarts && $scds[$i]->[3] == $start);
	    } else {
	      $scds[$i]->[2] = 'Terminal' if ($numstops && $scds[$i]->[3] == $stop);
	    }
	    if (($i ==  ($numcds -1))&&($scds[$i]->[6] eq '+')) {
	      $scds[$i]->[2] = 'Terminal' if ($numstops && $scds[$i]->[4] == $stop);
	    } else {
	      $scds[$i]->[2] = 'First' if ($numstarts && $scds[$i]->[4] == $start);
	    }
	  }
	} 
      }
    } else {
      if ($addstop && @scds) {
	if ($tgff[6] eq "+") {
	  $scds[-1]->[4]+=3;
	} else {
	  $scds[0]->[3]-=3;
	}
      }
      if ($fits) {
	if ($numcds == 1) {
	  $scds[0]->[2] = 'Single';
	} else {
	  for (my $i = 0;$i<$numcds;$i++) {
	    $scds[$i]->[2] = 'Internal';
	    $scds[$i]->[2] = 'First' if $i == 0;
	    $scds[$i]->[2] = 'Terminal' if $i == ($numcds -1);
	  }

	}
      
      }
    }
    if ($removestop && @scds) {
      if ($tgff[6] eq "+") {
	$scds[-1]->[4]-=3;
      } else {
	$scds[0]->[3]+=3;
      }
    }
    ### ADD UTRS
    if($addutr){
      foreach my $ex (@sexons){
	my @exon=@$ex;
	my $tags=$exon[8];
	$tags=~s/\;?ID=[^;]+//;
	if ($exon[3]<$CDSstart && $exon[4]<$CDSstart){
	  push @{$transcript},[$exon[0],$exon[1],'UTR',$exon[3],$exon[4],$exon[5],$exon[6],$exon[7],$tags];
	}elsif($exon[3]>$CDSend && $exon[4]>$CDSend){
	  push @{$transcript},[$exon[0],$exon[1],'UTR',$exon[3],$exon[4],$exon[5],$exon[6],$exon[7],$tags];
	}elsif($exon[3]<$CDSstart && $exon[4]>=$CDSstart){
	  push @{$transcript},[$exon[0],$exon[1],'UTR',$exon[3],$CDSstart-1,$exon[5],$exon[6],$exon[7],$tags];
	}elsif($exon[3]<=$CDSend && $exon[4]>$CDSend){
	  push @{$transcript},[$exon[0],$exon[1],'UTR',$CDSend+1,$exon[4],$exon[5],$exon[6],$exon[7],$tags];
	}
      }
    }


    foreach my $gffrec (@{$transcript}) {
      if ($addname){
	if ($gffrec->[8] !~/Name=/){
	  if ($gffrec->[8] =~/ID=([^;]+)/){
	    $gffrec->[8].=";Name=$1";
	    $gffrec->[8]=~s/^;//;
	  }
	}
      }
      if ($ov != 2) {
	$gffrec->[8]=~s/;$//;
	print (join("\t",@$gffrec),"\n");
       } else {
       	if ($gffrec->[2] =~ m/exon|First|Internal|Terminal|Single|CDS|UTR/) {
       	  print (join("\t",@$gffrec),"\n");
       	}
       }
    }
    print "###\n" if ($ov!=2);
  }

print STDERR "done!\n";
sub gffsort
  {
    $a->[0] cmp $b->[0]
      ||
	$a->[3] <=> $b->[3]
	 ||
	 $a->[4] <=> $b->[4]
	  ||
	  $a->[6] <=> $b->[6]
}

sub cdsPhase {
  my($db, $cdsA)= @_;
  ## assume cdsA are all/only cds exon set for one gene/mrna
  return $cdsA unless(ref $db);
  
  my $cstrand= $cdsA->[0]->[6];
  my $isrev= ($cstrand eq '-' || $cstrand < 0);
  
  my @cds;   # sort by start
  if ($isrev) { @cds= sort{ $b->[3] <=> $a->[3] } @$cdsA; } # end 1st
  else { @cds= sort{ $a->[3] <=> $b->[3] } @$cdsA; } # start 1st
  my $nt_length= 0;
  my $ispartial= 0;
  
  foreach my $ix (0 .. $#cds)  {
    my($ref,$src,$type,$start,$stop,$score,$strand,$phase,$attr)= @{$cds[$ix]};
    my $id=$attr; $id =~ s/^(ID|Parent)=//; $id =~ s/[;,\s].*$//;
    # do we check exon ordering? use as given in gff?  need 1st .. last, differs for strands
    
    if($ix == 0) { 
      ## 1st exon; find start ATG; ** only need 3 bases at start, not all
      $ispartial= 0;
      my($bstart,$blen)= ($start - 1, $stop - $start + 1);
      #my($bstart,$blen)= ($isrev) ? ($stop-8,8) : ($start-1, 8);
      my $exondna = SeqOp::get_seq_BioDBFasta($db,$ref,$start,$stop,$strand);
      #my $exondna  = substr( $refdna->seq(), $bstart, $blen);
      #if($isrev) {  
      #  $exondna = reverse $exondna;
      #  $exondna =~ tr/gatcGATC/ctagCTAG/;
      #  }
      my $inc5= 0;
      for (; $inc5<=3; $inc5++) {
        my $atg= substr($exondna, $inc5, 3);
        last if($atg =~ /atg/i);
        }

      ## fixme, if $ispartial probably need check best aa translation frame
      ## yes; need full cds/all exons and translate() method
      $ispartial=1 if $inc5 || $phase;
      if ($ispartial) { 
        $nt_length = 0; $inc5 = 0; #start not found/incomplete prot ?
        my $cdsdna= $exondna;
        foreach my $ex (1 .. $#cds) {
          my($ref1,$src1,$type1,$start1,$stop1,$score1,$strand1,$phase1,$attr1)= @{$cds[$ex]};
	  my $exon2dna = SeqOp::get_seq_BioDBFasta($db,$ref1,$start1,$stop1,$strand1);
          #my($bstart,$blen)= ($start1 - 1, $stop1 - $start1 + 1);
          #my $exon2dna  = substr( $refdna->seq(), $bstart, $blen);
          #if($strand1 eq '-' || $strand1 < 0) {  
          #  $exon2dna = reverse $exon2dna;
          #  $exon2dna =~ tr/gatcGATC/ctagCTAG/;
          #  }
          $cdsdna.= $exon2dna;
          }
        $inc5 = getBestFrame( $cdsdna, $id, $phase);
        }
      
      if ($inc5 == 1) { $nt_length = 2; }
      elsif ($inc5 == 2) { $nt_length = 1; }
      else  { $nt_length = 0; }
    }
    
    my($inc5,$inc3,$elength,$frame);
    $elength = $stop - $start + 1;
		$nt_length  += $elength;
		$inc3        = $nt_length % 3;
		$inc5        = ($elength - $inc3) % 3; # only care about this one
		$frame       = ($start + $inc5) % 3;
		if ($inc5 == -1) { $inc5 = 2; }
    
    my $changed=0;
    if ($phase eq '.') {  $changed=1; }
    elsif ($phase ne $inc5 ) { 
      $changed=2; 
      warn "# phase change exon[$ix]: $phase => $inc5; $ref:$start-$stop/$strand,$type:$src,$id\n" if $verbose;
      } 
    if($changed) { $cds[$ix]->[7]= $inc5; }  
    if($ispartial && $ix == 0) { $cds[$ix]->[8] .= ";partial_gene=true"; } # 5prime_partial=true; 3prime..
    }
    
  return \@cds;
}



my @s5CodonTable = ();
BEGIN{
 @s5CodonTable = (
	 [
		 ['K','N','K','N','X',],
		 ['T','T','T','T','T',],
		 ['R','S','R','S','X',],
		 ['I','I','M','I','X',],
		 ['X','X','X','X','X',],
	],
	 [
		 ['Q','H','Q','H','X',],
		 ['P','P','P','P','P',],
		 ['R','R','R','R','R',],
		 ['L','L','L','L','L',],
		 ['X','X','X','X','X',],
	],
	 [
		 ['E','D','E','D','X',],
		 ['A','A','A','A','A',],
		 ['G','G','G','G','G',],
		 ['V','V','V','V','V',],
		 ['X','X','X','X','X',],
	],
	 [
		 ['*','Y','*','Y','X',],
		 ['S','S','S','S','S',],
		 ['*','C','W','C','X',],
		 ['L','F','L','F','X',],
		 ['X','X','X','X','X',],
	],
	 [
		 ['X','X','X','X','X',],
		 ['X','X','X','X','X',],
		 ['X','X','X','X','X',],
		 ['X','X','X','X','X',],
		 ['X','X','X','X','X',],
	],

);
}

sub ibase {
  my $c= substr($_[0],$_[1],1);
  return 0 if ($c eq 'A');
  return 1 if ($c eq 'C');
  return 2 if ($c eq 'G');
  return 3 if ($c eq 'T');
  return 4;
}  
  
sub translate {
  my($cds, $offset)= @_;
  $cds = uc($cds); ## fix chars ??
  my $aa="";
  my $aa_length = int((length($cds) - $offset) / 3);
	for (my $i = 0; $i < $aa_length; $i++) {
		my $idx = 3 * $i + $offset;
		$aa .= $s5CodonTable[ ibase($cds,$idx)][ ibase($cds,$idx+1) ][ ibase($cds,$idx+2) ];
	}
  return $aa; 
}

sub getBestFrame {
  my($cds, $id, $original_phase)= @_;
  my $orig_score = -999;
  my ($bestscore,$besti)= (-999,0);
  for (my $i= 0; $i<3; $i++) {
    my $pro= translate( $cds,$i );
    $cds=~s/(TAA|TAG|TGA)$//i;
    my $coding_potential=compute_cp($cds,(3-$i)%3);
    #print STDERR "$i\t$coding_potential\n";
    my $score = $pro =~ tr/*/*/; # has_internal_stops($pro);
    #is inner M bad?# $score += $pro =~ tr/M/M/;   # has_internal_starts($pro);
    $score *= -3; 
    if (substr($pro,length($pro)-1,1) eq '*') { $score += 4; } # adj internal == end
    #if (substr($pro,0,1) eq 'M') { $score += 1; }
    $score+=$coding_potential;
    if ($i == $original_phase){$orig_score = $score;}
    if ($score >= $bestscore) { 
      if($score == $orig_score){ 
	$besti= $original_phase; 
      }else{
	$besti= $i;
      }
      $bestscore=$score;
    }
 warn("# bestFrame[$i,$id]: $score ; $pro \n") if($verbose); # debug  
    }
  return $besti;
}

sub compute_cp {
  my $seq   = shift;
  my $ucseq = uc($seq);
  #print STDERR $ucseq,"\n";
  my $frame = shift;
  my $f = $frame;
  my $score = 0;
  for(my $p=0;$p<(length($ucseq)-5);$p++){
    my $h=substr($ucseq,$p,6);
    #print STDERR "$h\t$f\t$llhex{$h}{$f}\n";
    $score+=$llhex{$h}{$f};
    $f++;
    if($f==3){$f=0;}
  }
  return $score;
}

