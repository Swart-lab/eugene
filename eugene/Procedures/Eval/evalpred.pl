#!/usr/bin/perl

# Programme qui evalue la precision d'une serie 
# d'annotation de sequences nucleiques d'EuGene,
# prend en compte les predictions reverse,
# accepte plusieurs genes par sequence soumise, ou par prediction.

# arguments :
# 1 : coordonnees reelles des exons
# -un espace avant chaque coordonnee exon
# -un gene par ligne
# -ligne vide entre chaque sequence soumise
# -et deux lignes vides en fin du fichier!!
# exemple format coord:

#  125 653 2789 3052
#  3238 3267 3289 3302
#
#  13 56 1875 1945 6211 6574
#
#  26 38 348 392 ...
#  ...
#
# (et possible -23 -48 -139 -156...)

# 2 : coordonnees predites.
# (sorties EuGeneAS);(ASclean, algogen_coord.out, ..)
# exemple format predictions:

#  Init  +   1   13   13   +1   +1   0     14   1.0
#  Term  +  267 1282 1016  +2   +2   266   1283  1.
#  Init  + 1394 1437 ...
#  Intr ...
#  Term ...
#
#  Init  + 1111 1120  88   +1   +3  1111  1120 1.0
#  Sngl  - 1890 1945  ...
# ...
# depuis mai 2002: tolere les lignes "Utr5" et "Utr3"

# 3 facultatif (immediatement precede de -o)
# offset, qui definit des bornes de part et d'autre des 
# coordonnees des genes reels au dela desquelles les predictions
# ne seront pas prises en compte. (cf evaluation avec araset, ou offset=300)

$setoffset=0;
$offset=0;

if ($#ARGV!=1){
  ($ARGV[2]=~ s/\-o//) || die "Usage : <fichier coordonnees reelles des exons> <fichier coordonnees predites> (et eventuellement<-oN>, N=offset)\n";
  $setoffset=1;  
  $offset= (abs($ARGV[2]));
  print "offset= $offset\n";
}

### Fonction d'evaluation au niveau gene qui,
# pour une sequence donnee, regarde si chaque
# gene reel a une prediction parfaitement egale
sub evaluation_niveau_gene{
  my@COORD2D=@{@_[0]}; # 1er argument
  my@PRED2D=@{@_[1]};  # 2eme argument

  # nbre de genes reels/predits
  $nGR+=($#COORD2D + 1);
  $nGP+=($#PRED2D + 1);

  # Vrais positifs:
  # pour chaque gene reel de la sequence
  foreach $pointeur_gene_reel (@COORD2D){
    my@genereelcoord = @{$pointeur_gene_reel};
    foreach $pointeur_gene_pred (@PRED2D){
      if("@genereelcoord" eq "@{$pointeur_gene_pred}"){
	$VPg++;
      }
    }
  }
}

## Fonction d'evaluation au niveau exon
sub evaluation_niveau_exon{
  my@E=@{@_[0]};  # 1er argument
  my@P=@{@_[1]};  # 2eme argument
  my$i=0;;
  my$j=0;;

  # nbre d'exons reels et predits:
  $nER+=($#E+1)/2;
  $nEP+=($#P+1)/2;

  # Vrais positifs:
  for($i=1;$i<=$#E;$i+=2){ # toutes les fins d'exons reels
    for($j=1;$j<=$#P;$j+=2){ # toutes les fins d'exons predits
      # exon correctement predit:
      $VPe+=(($E[$i-1]==$P[$j-1])&&($E[$i]==$P[$j]));
    }
  }
}

########    "LEGENDE" ou guide des symboles ###################
##
##  FP/FN=Faux positifs/negatifs,
##  VP/VN=Vrais positifs/negatifs,
##  suivi de g=gene, e=exon, ou n=nucleotide.
##
##  tableau E =coordonnees des exons reels pour le gene courant
##  tableau P =coordonnees des exons predits pour le gene courant
##
##  nG/nE/nN = nombre de genes, exons et nucleotides,
##  suivi de R=reellement codants, ou P=predits comme codants.
##  
################################################################

$FPg=0;$FNg=0;$VNg=0;$VPg=0;$nGR=0;$nGP=0;
$FPe=0;$FNe=0;$VNe=0;$VPe=0;$nER=0;$nEP=0;
$FPn=0;$FNn=0;$VNn=0;$VPn=0;$nNR=0;$nNP=0;$nNtot=0;
@coord=();@pred=();
@E=();
@P=();
@COORD2D=();
@PRED2D=();
$nseqreelles=1;
$nseqpred=1;

open(COORD,"$ARGV[0]") || die "Can't open fich coord reelles $ARGV[0]\n";
open(PRED,"$ARGV[1]") || die "Can't open fich coord predites $ARGV[1]\n";
#open(LG,"$ARGV[2]") || die "Can't open liste des longueurs seq $ARGV[2]\n";

$nlCOORD=0;$nlPRED=0;
$i=0;$j=0;


#print("SEQUENCE REELE n°$nseqreelles\n");
while($lCOORD=<COORD>){
  chomp $lCOORD;
  $nlCOORD++;

  # Si on est dans une sequence
  if($lCOORD =~ /(\s+\-?([0-9]+))+/){
    $i++;
    @TMP=split(/\s+/,$lCOORD);
    shift @TMP;  # case vide
    foreach $coord (@TMP){
      push (@E,$coord); #stocke ensemble des exons de la sequence
    }
    @COORD2D[$i-1]=([@TMP]);   #stocke ensemble des genes de la sequence
  }
  else{
    if($lCOORD == ""){ # Fin de la sequence

      # calcul des bornes droite et gauche
      $borneg= abs($COORD2D[0][0]) - $offset;
      $borned= abs($COORD2D[$#COORD2D][$#{$#COORD2D}]) + $offset;

      ## LECTURE des predictions d'EuGene
      $j=0;
      $gene_en_cours=0;
      foreach $pointeur2D (@PRED2D){
	@{$pointeur2D}=(); }
      @PRED2D=();
      @P=();
      @TMP=();

      chomp($lPRED=<PRED>);$nlPRED++;
      if($lPRED eq ""){ # pas de gene trouve
	$j++;
#	chomp($lPRED=<PRED>); # direction prochaine ligne vide de separation
      }

      while($lPRED ne ""){ # tant qu'il y a un exon trouve
	if ($lPRED =~ /\s+([a-zA-Z]+)\s+([\+\-])\s+([0-9]+)\s+([0-9]+)\s+[0-9]+\s+[^\s]+\s+[^\s]+\s+([0-9]+)\s+([0-9]+)\s+/) {
	  if ( ($setoffset==0) ||
	       (($4 >= $borneg) && ($3 <= $borned))) {
	    $type=$1;
	    $sens=$2;
	    $gene_en_cours=1;
	    # Si format ou fichier mauvais:
	    if (!( ($type eq "Init")||($type eq "Intr")||($type eq "Term")||($type eq "Sngl"))){ die" Pb format fichier predictions $ARGV[1] ligne $nlPRED\n"};
	    
	    # Stockage des exons
	    $signe=( ("$sens" eq '+')? "" : '-');
	    push(@P,"$signe$3");
	    push(@P,"$signe$4");
	    push(@TMP,"$signe$3");
	    push(@TMP,"$signe$4");
	    
	    ## fin d'un gene
	    if ( (($type eq "Term")&&("$sens" eq '+'))||(($type eq "Init")&&("$sens" eq '-'))||($type eq "Sngl")){
	      $j++;
	      $gene_en_cours=0;
	      @PRED2D[$j-1]=[@TMP];
	      @TMP=();
	    }
	  }
	}
	elsif($lPRED !~ /\s+Utr\d\s/) {die"Pb format fichier predictions $ARGV[1] ligne $nlPRED\n"};
	chomp($lPRED=<PRED>);
	$nlPRED++;
      }
      if ($gene_en_cours==1) {
	# Cas ou la derniere prediction "continue" apres la sequence
	# (gene non termine par un "+Term" ou "-Init")
	$j++;
	@PRED2D[$j-1]=[@TMP];
	@TMP=();
	$gene_en_cours=0;
      }

      &evaluation_niveau_gene(\@COORD2D,\@PRED2D);
      &evaluation_niveau_exon(\@E,\@P);
      $i=0;

      # liberation memoire tableaux:
      @E=();
      @P=();
      foreach $pointeur2D (@COORD2D){
	@{$pointeur2D}=(); }
      @COORD2D=();
      foreach $pointeur2D (@PRED2D){
	@{$pointeur2D}=(); }
      @PRED2D=();
    }
    else{
      die"Pb format fichier coordonnees $ARGV[0] ligne $nlCOORD\n";
    }
  }
}

if (($nGR==0)||($nER==0)){die"Pb : nbre de seq soumises (ou d'exons) nul!\n"}

$SENSg=($VPg*100)/$nGR;
$SENSe=($VPe*100)/$nER;
#$SENSn=($VPn*100)/$nNR;

if (($nGP==0)||($nEP==0)){
  $SPECg=100;
  $SPECe=100;
#  $SPECn=100;
}
else{
  $SPECg=($VPg*100)/$nGP;
  $SPECe=($VPe*100)/$nEP;
#  $SPECn=($VPn*100)/$nNP;
}

print "\n>TOTAL\n";
print "$VPg genes bien detectes sur $nGR avec $nGP predictions\n";
print "$VPe exons bien detectes sur $nER avec $nEP predictions\n";
#print "$VPn nt codants bien detectes sur $nNR ($FNn rates);$nNP predits dont $FPn faux pos\n";
#print "\n";
#print "SENSIBILITE GENES : $SENSg\n";
#print "SPECIFICITE GENES : $SPECg\n";
#print "SENSIBILITE EXONS : $SENSe\n";
#print "SPECIFICITE EXONS : $SPECe\n";
#print "SENSIBILITE NUCL  : $SENSn\n";
#print "SPECIFICITE NUCL  : $SPECn\n";

#print"$SENSg $SPECg $SENSe $SPECe $SENSn $SPECn\n";
print"SNG: $SENSg SPG: $SPECg SNE: $SENSe SPE: $SPECe\n";

close COORD;
close PRED;
