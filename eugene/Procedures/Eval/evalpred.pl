#!/usr/bin/perl

$setoffset=0;
$offset=0;
$sortie=1; # 0-> short, 2-> detaillee
$name2="";
$NOMduPROG = $0;
@CMD="$NOMduPROG "."@ARGV";

sub description{
  die <<EOF

############################################################################ 
# Ce prog evalue la qualite de predictions d'EuGene sur des seq nucleiques,
# gere les predictions reverses et plsrs genes ou predictions par sequence. 
#
# arguments :
#
#     - 1 : fichier coordonnees reelles des exons
# -les sequences annotees sont dans le meme ordre que les predictions
# -un espace devant chaque coordonnee exon
# -un signe moins devant la coord si le sens est reverse
# -un gene par ligne
#  (retour chariot apres le dernier exon d'un gene)
# -une ligne vide apres chaque sequence (derniere comprise)
#  (ligne vide = seulement retour chariot "\n")
# exemple format coord:
#
#  13 56 1875 1945 6211 6574
#
#  125 653 2789 3052
#  -3238 -3267 -3289 -3302
#
#  76 331
#  ...
#
#     - 2 : coordonnees predites.
# sortie standard d'EuGeneAS 
# ex: EuGeneAS *.fasta > out, ou EuGeneAS `cat seqfiles.list` > out
# exemple format predictions:
#
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
#
#     - 3, 4 et 5: arguments facultatifs (ordre indifferent).
# -oN
# N=offset, qui definit des bornes de part et d'autre des 
# coordonnees des genes reels au dela desquelles les predictions
# ne seront pas prises en compte.(cf evaluation avec araset, ou offset=300)
# -pX
# affichage de la sortie, X=s pour short (utilise pour l optimisation),
# et X=l pour long (equivaut a sortie detaillee) (par defaut:intermediaire)
# -fX
# sert pour sortie detaillee, le nom des sequences proviendra non pas du
# fichier predictions (comme par defaut, preleve sur la 1ere colonne),
# mais du fichier X (un nom par ligne, pouvant etre l ID, le chemin,...)
#
# ex : $NOMduPROG araclean.exons.coord eugene.outpred -o300 -pl -fID.list
#
# todo: ajouter niveau nucleotidique 
# (demande un autre fichier en entree avec les lgrs des seq)
############################################################################

EOF
}

############################################################################ 
sub usage{
  die <<EOF

usage : 
$NOMduPROG exons_coord_file pred_file [-o{offset}] [-p{s|l}] [-f{IDfile}]
        -o  : offset, only predictions in the region 
              "gene +/- offset" are considered.
        -ps : short outpout (used for optimisation scripts)
        -pl : long outpout  (used for expertised analyses)
              (default = intermediaire)
        -f  : name (or ID, path...) file (one name/line), used if -pl active 
"-" arguments are optionnal, must be adjacent to the "-", order no matters
ex : $NOMduPROG seq.exons.coord eugene.outpred -o300 -pl -fseq.ID.list

EOF
}

############################################################################ 
# Fonction d'evaluation au niveau gene qui,
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

############################################################################ 
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

############################################################################ 
# Renvoie le nombre d'exons identiques entre 2 genes
sub nombre_exons_identiques{
  my@gene1=@{@_[0]};
  my@gene2=@{@_[1]};
  my$i=0;
  my$j=0;
  my$neID=0;
  for($i=1;$i<=$#gene1;$i+=2){ # toutes les fins d'exons gene1
    for($j=1;$j<=$#gene2;$j+=2){ # toutes les fins d'exons gene2
      # exon identique:
      $neID+=(($gene1[$i-1]==$gene2[$j-1])&&($gene1[$i]==$gene2[$j]));
    }
  }
  return $neID;
}

############################################################################ 
# Fonction nouvelle d'evaluation tous niveaux
# et d'affichage des resultats
sub evaluation_et_affichage{
  my@COORD=@{@_[0]}; # 1er argument
  my@PRED=@{@_[1]};  # 2eme argument
  my@E=@{@_[2]};
  my@P=@{@_[3]};
  my$seq=@_[4];

  # nbre de genes reels/predits
  my$ngrtot=($#COORD + 1);
  my$ngptot=($#PRED + 1);
  $nGR+=$ngrtot;
  $nGP+=$ngptot;

  # nbre d'exons reels et predits:
  my$nertot=($#E+1)/2;
  my$neptot=($#P+1)/2;
  $nEP+=$neptot;
  $nER+=$nertot;

  my$i=0;
  my$j=0;
  my$c="";
  my$p="";
  my$npredparfaites=0;
  my$nexonstrouves=0;
  print("\n -SEQUENCE ANALYSEE: $seq\n");

  my$ngr=0;
  my$ngp=0;
  @PRINT=();

  # pour chaque gene reel de la sequence
  foreach $c (@COORD){
    $ngr++;
    my@genereel = @{$c};
    my$ne=($#genereel+1)/2;
    # affichage coordonnees reelles du gene:
    print("gene reel  $ngr/$ngrtot:");
    for($i=0;$i<$#genereel;$i+=2){
      print(" $genereel[$i]-$genereel[$i+1]");
    }
    print(' ... ');
    $predit=0;
    $ngp=0;
    foreach $p (@PRED){
      # on cherche parmi les genes predits
      $ngp++;
      if("@genereel" eq "@{$p}"){
	# Vrais positifs genes:
	$VPg++;
	$npredparfaites++;
	$predit=2;
	print("OK! (detecte en pred n°$ngp)");
      }
      else{
	$neid=&nombre_exons_identiques(\@genereel,$p);
	if ($neid>0){
	  print("rate! ($neid exons trouve(s) sur $ne en pred n°$ngp)");
	  $predit=1;
	}
      }
    }
    ($predit<2) ? push(@PRINT,"fng") : push(@PRINT,"vpg");
    ($predit>=1)? print("\n") : print("tous les exons rates!\n");
  }

  # affichage coordonnees des genes predits:
  $ngp=0;
  foreach $p (@PRED){
    $ngp++;
    my@genepred=@{$p};
    my$ne=($#genepred+1)/2;
    print("prediction $ngp/$ngptot:");
    for($i=0;$i<$#genepred;$i+=2){
      print(" $genepred[$i]-$genepred[$i+1]");
    }
    print(" ... ");
    $predok=0;
    $numerogene=0;
    foreach $p (@COORD){
      $numerogene++;
      if("@genepred" eq "@{$p}"){
	$predok=2;
	print("OK! (parfaite pour gene n°$numerogene)");
      }
      else{
	$neid=&nombre_exons_identiques(\@genepred,$p);
	if ($neid>0){
	  print("imparfaite (cf gene n°$numerogene)");
	  $predok=1;
	}
      }
    }
    if ($predok<2) { push(@PRINT,"fpg") } # vpg deja mis
    ($predok==0)? print("faux positif!\n") : print("\n");
  }

  # EXONS:
  $tmpe=&nombre_exons_identiques(\@E,\@P);
  $VPe+= $tmpe;
  $nexonstrouves+= $tmpe;

  for($i=1;$i<=$#P;$i+=2){ # toutes les fins d'exons P
    my$flag=0;
    for($j=1;$j<=$#E;$j+=2){ # toutes les fins d'exons E
      # exon identique:
      $flag+=(($E[$i-1]==$P[$j-1])&&($E[$i]==$P[$j]));
    }
    ($flag==0) ? push(@PRINT,"fpe") : push(@PRINT,"vpe");
  }
  for($i=1;$i<=$#E;$i+=2){ # toutes les fins d'exons E
    my$flag=0;
    for($j=1;$j<=$#P;$j+=2){ # toutes les fins d'exons P
      # exon identique:
      $flag+=(($E[$i-1]==$P[$j-1])&&($E[$i]==$P[$j]));
    }
    if ($flag==0) { push(@PRINT,"fne") }
  }

  print("Code  pour $seq: @PRINT\n");

  print("Total pour $seq: ");
#  print("Pour $ngrtot genes reels, $ngptot prediction(s), dont $npredparfaites parfaite(s) \n");
#  print("Pour $nertot exons reels, $neptot prediction(s), dont $nexonstrouves parfaite(s)\n\n");
#  print("$npredparfaites genes bien predits sur $ngrtot (avec $ngptot predictions)\n");
#  print("$nexonstrouves exons bien predits sur $nertot (avec $neptot predictions)\n");
  print("$ngrtot gene(s), $ngptot predit(s), $npredparfaites trouve(s), $nertot exon(s), $neptot predit(s), $nexonstrouves trouve(s)\n");
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

######################################################################
# lecture des arguments (pas joli, mais pas besoin de use Getopt::Std)
if ($#ARGV != 1){
  if ($ARGV[0] eq '-h') { description() }
  if ($#ARGV < 1) { usage() }
  for (my$i=2 ; $i<=$#ARGV ; $i++) {
    if ($ARGV[$i] =~ /^\-(.)(.+)/) {
      if ($1 eq 'o') {
	if (!($2 =~ /^\d+$/)) { usage() }
	$setoffset=1;
	$offset= $2;
      }
      elsif ($1 eq 'f') {
	$name2="1";
	open(LS,"$2") || die "Can't open seqID list $2\n";
      }
      elsif ($1 eq 'p') {
	if (($2 eq 's') || ($2 eq 'l')) {
	  $sortie= ( ($2 eq 's') ? 0 : 2 );
	}
	else { usage() }
      }
      else { usage() }
    }
    else { usage() }
  }
}

open(COORD,"$ARGV[0]") || die "Can't open fich coord reelles $ARGV[0]\n";
open(PRED,"$ARGV[1]") || die "Can't open fich coord predites $ARGV[1]\n";

$nlCOORD=0;$nlPRED=0;
$numgene=0;$j=0;
$newseqreelle=1;

if ($sortie!=0) {
  print("\nEVALUATION DES PREDICTIONS D'EUGENE (@CMD)\n");
}

while($lCOORD=<COORD>) {
  chomp $lCOORD;
  $nlCOORD++;
  if (($name2) && ($newseqreelle==1)) {
    chomp($name2=<LS>);
    $newseqreelle=0;
  }
  # Si on est dans une sequence
  if($lCOORD =~ /(\s+\-?([0-9]+))+/){
    $numgene++;
    @TMP=split(/\s+/,$lCOORD);
#TMP!! on pourrait verifier si vide avant de shifter (+ de souplesse de format coord)
    shift @TMP;  # case vide
    foreach $coord (@TMP){
      push (@E,$coord); #stocke ensemble des exons de la sequence
    }
    @COORD2D[$numgene-1]=([@TMP]);   #stocke ensemble des genes de la sequence
  }
  else{
    if($lCOORD == ""){ # Fin de la sequence
      $newseqreelle=1;

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
	$name= ( ($name2) ? $name2 : " -NO_NAME- " );
      }

      while($lPRED ne ""){ # tant qu'il y a un exon trouve
	if ($lPRED =~ /^([^\s]+)\s+([a-zA-Z]+)\s+([\+\-])\s+([0-9]+)\s+([0-9]+)\s+[0-9]+\s+[^\s]+\s+[^\s]+\s+([0-9]+)\s+([0-9]+)\s+/) {
	  if ( ($setoffset==0) ||
	       (($5 >= $borneg) && ($4 <= $borned))) {
	    $name= ( ($name2) ? $name2 : $1);
	    $type=$2;
	    $sens=$3;
	    $beg=$4;
	    $end=$5;
#	    ($name= $1) =~ s/\.\d+\.\d+\.\d+//;
	    $gene_en_cours=1;
	    # Si format ou fichier mauvais:
	    if (!( ($type eq "Init")||($type eq "Intr")||($type eq "Term")||($type eq "Sngl"))){ die" Pb(1) format fichier predictions $ARGV[1] ligne $nlPRED\n"};

	    # Stockage des exons
	    $signe=( ("$sens" eq '+')? "" : '-');
	    push(@P,"$signe$beg");
	    push(@P,"$signe$end");
	    push(@TMP,"$signe$beg");
	    push(@TMP,"$signe$end");

	    ## fin d'un gene
	    if ( (($type eq "Term")&&("$sens" eq '+'))||(($type eq "Init")&&("$sens" eq '-'))||($type eq "Sngl")){
	      $j++;
	      $gene_en_cours=0;
	      @PRED2D[$j-1]=[@TMP];
	      @TMP=();
	    }
	  }
	}
	elsif($lPRED !~ /\s+Utr\d\s/) {die"Pb(2) format fichier predictions $ARGV[1] ligne $nlPRED\n"};
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
      if ($sortie == 2) {
	&evaluation_et_affichage(\@COORD2D,\@PRED2D,\@E,\@P,$name);
      }
      else {
	&evaluation_niveau_gene(\@COORD2D,\@PRED2D);
	&evaluation_niveau_exon(\@E,\@P);
      }
      $numgene=0;

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
      die"Pb(3) format fichier coordonnees $ARGV[0] ligne $nlCOORD\n";
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
if ($sortie==0) {
  print"$SENSg $SPECg $SENSe $SPECe\n";
}
else {
  print "\n>TOTAL (@CMD)\n";
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
  print"SNG: $SENSg SPG: $SPECg SNE: $SENSe SPE: $SPECe\n";
}
if($name2) {close LS}
close COORD;
close PRED;
