#!/usr/bin/perl

$usage .= "$0 -- convert coord format to General Feature Format (GFF)\n";
$usage .= "Infos about coord format : see evalpred.pl script\n";
$usage .= "Infos about GFF format  : see
 http://www.sanger.ac.uk/Software/formats/GFF/GFF_Spec.shtml\n";
$usage .= "Cannot handle multi-sequences coord files (one seq per .coord file)\n";
$usage .= "Can handle several genes in the seq\n";
$usage .= "doesn't compute the frame\n";
$usage .= "\n";
$usage .= "Usage: $0 [input coord file]\n";
$usage .= "\n";

$Init="E.Init";
$Intr="E.Intr";
$Term="E.Term";
$Single="E.Single";

if ($#ARGV != 0) {
  die "$usage";
}

$coordfile = "$ARGV[0]";

open(COORD,"$coordfile") || die "Can't open input coord file $coordfile\n";

$ngene=0;
while (<COORD>){
  if (/\d/) {
    $ngene++;
    @coord=split;
    if ($coord[0] =~ /^\-?\d+$/) { # first column = coordinate, no seq name (old coord format)
      $name="$coordfile".".$ngene";
    }
    else {
      $name= "$coord[0]".".$ngene";
      shift(@coord);
    }

    $strand= ( ($coord[0] <0) ? '-' : '+' );
    $coord =~ s/-//;
    foreach (@coord) { s/-// }
    for ($i=0;$i<$#coord;$i+=2) {
      if ($#coord==1) {
	$feature="$Single";
      }
      elsif ( ($i==0) || ($i==$#coord-1)) {
	if ($i==0) {
	  $feature= ( ($strand eq '+') ? "$Init" : "$Term" );
	}
	else {
	  $feature= ( ($strand eq '+') ? "$Term" : "$Init" );
	}
      }
      else {
	$feature= "$Intr";
      }
      print "$name\t$0\t$feature\t$coord[$i]\t$coord[$i+1]\t\t.\t$strand\t.\n";
    }
  }
}

