#!/usr/bin/perl
# from a single multi-sequences coord file (e.g. araclean.coord)
# with the sequences list file (e.g. araclean.fasta)
# generates one .coord file per sequence

# WARNING: ERASE preceding .coord file if $ext == ".coord"!
# coord file must contain only complete genes

($#ARGV==1) || die "Usage: $0 coordfile fastalistfile;\nWARNING: ERASE preceding .coord files!!\n          coord file must contain only complete genes\n";

$ext='.coord';
$nseq=0;

$coordfile=@ARGV[0];
$list=@ARGV[1];

open(COORD,"$coordfile")||die"can't open coord file $coordfile\n";
open(SEQ,"$list")||die"can't open list file $list\n";

while($seq=<SEQ>) {
  chomp $seq;
  ($ID = $seq) =~ s/.+\///g;     # all after the last "/"
  print "\n-sequence $seq\n";
  $n=0;
  $coordfilename= "$seq"."$ext";
  open(OUT,">$coordfilename")||die"can't open coordfile $coordfilename\n";
  while($lcoord=<COORD>) {
    if ($lcoord =~ /(\s+\-?\d+)+/) {
      $n++;
      print OUT "$lcoord";
    }
    elsif ($lcoord == "") {
      last;
    }
  }
  close OUT;
}
close COORD;
close SEQ;
