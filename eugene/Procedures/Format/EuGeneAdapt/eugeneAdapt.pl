#!/usr/bin/perl

use strict;
use warnings;
use IO::Handle;
use Getopt::Long;
use gffUtils;

#=============================================================================#
#=             Copyright (c) 2003 by INRA. All rights reserved.              =#
#=                 Redistribution is not permitted without                   =#
#=                 the express written permission of INRA.                   =#
#=                     Mail : tschiex@toulouse.inra.fr                       =#
#=---------------------------------------------------------------------------=#
#= File         : eugeneAdapt.pl                                             =#
#= Description  : First step for the adaptation of EuGene to a new specie.   =#
#= Authors      : P.Bardou, S.Foissac, M.J.Cros, A.Moisan, T.Schiex          =#
#= History      : version 1.0 (Dec. 1, 2002)  -> sim2eugene.pl   	     =#
#=                version 2.0 (Aou. 5, 2003)  -> eugeneAdapt.pl              =#
#= Improvements :                                                            =#
#=============================================================================#
# Warning :
#  - This script doesn't accept redondances in each list. Thus, if one genomic
#    sequence "contains" more than one cDNA (2), or if the genomic and the gff
#    file contains more than one gene, it has to be splited.
#    The alphabetical order allow to lauch EuGene using *.fasta and obtain the
#    predictions must be in the same order than the list files.
#  - For the second approch, each cDNA has to be the real cDNA corresponding
#    to the gene present in his genomic sequence, a good sequence quality is
#    request (> ~95%), in order to obtain with the sim4 alignment the true
#    exons coordinates.
#  - eugeneAdapt need :
#       -> two EMBOSS tools : extractseq and revseq.
#       -> Eu-imm2UTR       : Matrices builder
#       -> WAMbuilder       : WAM builder
#       -> seqlogo          : Plot signal consensus (WAM)
#       -> For the second approch sim4.
#  - eugeneAdapt creates temporary files deleted at the end of procedure.
#=============================================================================#

=pod

=head1

=head1 DESCRIPTION

  eugeneAdapt creates the different files need for an adaptation of
  EuGene to a new specie.
  Two approaches could be used. The first one (1) is based  on  the
  gene  coordinates  knowledge about  genomic  sequences.  And  the
  second (2) one is based on sim4  program (cDNA  and  Genomic  DNA
  Alignment).

=head1 OPTIONS

  -g, --genomic [FILE NAME] (1) & (2)
       File containing a list of the fasta genomic sequence file(s)
       sorted by alphabetical order.

  -c, --coord [FILE NAME] (1)

  -m, --mRna [FILE NAME] (2)
       File containing a list of the  fasta  cDNA  file(s), a  cDNA
       from a line corresponding to the  DNA  from  the  same  line
       number in the file genomic list.

  -s, --ssCoord [FILE NAME] (2)
       File containing  ATG  and  STOP for each cDNA,  in  the same
       order than the others lists (separated by a tabulation).

  -o, --offset [INT] (1) & (2)
       Optional: offset (default = 1000) max number  of  nucleotids
       on each side of the START and STOP position.

  -w, --wam [INT] (1) & (2)
       Optional: build  the  WAM  files  and  test  the  new  fasta
       sequences.

  -d, --dirOutput [DIR NAME] (1) & (2)
       Optional: output directory (default = EuAdaptOutput).

=head1 EXAMPLES

  (1) eugeneAdapt -g genoListFile -c coordFile -w
  (2) eugeneAdapt -g genoListFile -m cDnaListFile -s ATGStopFile -w

=head1 Copyright (c) 2003 by INRA. All rights reserved.

=head1 Mail : tschiex@toulouse.inra.fr

=head1

=cut


# Program used
my $cmd_sim4       = "sim4";       # cDNA and Genomic DNA Alignment
my $cmd_extractseq = "extractseq"; # EMBOSS
my $cmd_revseq     = "revseq";     # EMBOSS
# Matrices builder
my $cmd_IMM        = "../../../SensorPlugins/MarkovIMM/GetData/TrainIMM";
# WAM builder
my $cmd_WAMBuilder = "../../../SensorPlugins/0_SensorTk/GetData/WAMbuilder";
my $cmd_seqlogo    = "seqlogo";    # seqlogo

# For the ouput files (extension or file name)
my $sim4_out_ext = ".sim4.A3";
my $newdna_ext   = ".fasta";
my $gff_ext      = ".fasta.gff";
my $coordFile    = "coord.txt";
my $exonFile     = "exon.txt";
my $intFile      = "intron.txt";
my $utr5File     = "utr5.txt";
my $utr3File     = "utr3.txt";
my $matFile      = "matrice.mat";
my $matLogFile   = "matrice.log";
my $wam_sta_fp   = "WAM.START.FP";
my $wam_don_fp   = "WAM.DON.FP";
my $wam_acc_fp   = "WAM.ACC.FP";
my $wam_sto_fp   = "WAM.STOP.FP";
my $wam_sta_tp   = "WAM.START.TP";
my $wam_don_tp   = "WAM.DON.TP";
my $wam_acc_tp   = "WAM.ACC.TP";
my $wam_sto_tp   = "WAM.STOP.TP";
my $verif        = "verif.txt";

# For print info
my $excluSeq      = 7;       # first method (seq excluded when exon < N nt)
my $nbExclu       = 0;       # number of excluded sequence
my $nbNoUTR5      = 0;       # number of sequence without utr5
my $nbNoUTR3      = 0;       # number of sequence without utr3
my $nbNoIntron    = 0;       # number of sequence without intron
my $nbSeq         = 0;       # number of sequence
my $nbSeqRev      = 0;       # number of reverse sequence
my $exonMinLen    = 100000;  # min length for exon
my $exonMaxLen    = 0;       # max length for exon
my $intronMinLen  = 100000;  # min length for intron
my $intronMaxLen  = 0;       # max length for intron

# For WAM
my $worder  = 1;    # WAM order
my $wamstaB = 3;    # Context lenght, B before signal
my $wamstaA = 3;    # Context lenght, A after  signal
my $wamdonB = 3;    # Context lenght, B before signal
my $wamdonA = 4;    # Context lenght, A after  signal
my $wamaccB = 3;    # Context lenght, B before signal
my $wamaccA = 3;    # Context lenght, A after  signal
my $wamstoB = 8;    # Context lenght, B before signal
my $wamstoA = 8;    # Context lenght, A after  signal

# Options
my $genomic;
my $cDna;
my $ssCoord;
my $coord;
my $offset = 1000;
my $outputDir = "EuAdaptOutput";
my $wam;

MAIN:
{
  if(!GetOptions("genomic=s"   => \$genomic,
		 "mRna=s"      => \$cDna,
		 "ssCoord=s"   => \$ssCoord,
		 "coord=s"     => \$coord,
		 "offset=i"    => \$offset,
		 "wam"         => \$wam,
		 "dirOutput=s" => \$outputDir)) { exit(1); }

  ## Verify the options ##
  OptVerif();

  ## Open output files ##
  my $error = "no";
  open(COORD_OUT,">$coordFile") || die $error = $coordFile;
  open(EXON,  ">$exonFile")     || die $error = $exonFile;
  open(INTRON,">$intFile")      || die $error = $intFile;
  open(UTR5,  ">$utr5File")     || die $error = $utr5File;
  open(UTR3,  ">$utr3File")     || die $error = $utr3File;
  if(defined($wam)) {
    open(FPSTART, ">$wam_sta_fp") || die $error = $wam_sta_fp;
    open(FPDON,   ">$wam_don_fp") || die $error = $wam_don_fp;
    open(FPACC,   ">$wam_acc_fp") || die $error = $wam_acc_fp;
    open(FPSTOP,  ">$wam_sto_fp") || die $error = $wam_sto_fp;
    open(TPSTART, ">$wam_sta_tp") || die $error = $wam_sta_tp;
    open(TPDON,   ">$wam_don_tp") || die $error = $wam_don_tp;
    open(TPACC,   ">$wam_acc_tp") || die $error = $wam_acc_tp;
    open(TPSTOP,  ">$wam_sto_tp") || die $error = $wam_sto_tp;
    open(VERIF,   ">$verif")      || die $error = $verif;
  }
  if($error ne "no") { die "ERROR: Could not open file $error\n"; }

  ## Choose the approach
  if(defined($coord)) { GenoCoord(); }  # first  method
  else                { GenoSim4();  }  # second method

  ## Build matrices ##
  print"------------------------------------".
    "-------------------------------\n";
  print "Launch TrainImm (Build matrices)..";
  system "$cmd_IMM $exonFile $intFile -5 $utr5File -3 $utr3File > $matLogFile";
  system "mv $outputDir/exon.txt.bin $matFile";
  print "$matFile....done\n";

  ## Build WAM ##
  if(defined($wam)) {
    print"------------------------------------".
      "-------------------------------\n";
    print "Launch WAMBuilder (Build WAM)..";
    print "StartTF..";
    my $arg = "$worder ".($wamstaB+$wamstaA+3)." ".($worder+1);
    system "$cmd_WAMBuilder $arg $wam_sta_tp $wam_sta_tp. >& $wam_sta_tp.log";
    system "$cmd_WAMBuilder $arg $wam_sta_fp $wam_sta_fp. >& $wam_sta_fp.log";
    print "DonTF..";
    $arg = "$worder ".($wamdonB+$wamdonA+2)." ".($worder+1);
    system "$cmd_WAMBuilder $arg $wam_don_tp $wam_don_tp. >& $wam_don_tp.log";
    system "$cmd_WAMBuilder $arg $wam_don_fp $wam_don_fp. >& $wam_don_fp.log";
    print "AccTF..";
    $arg = "$worder ".($wamaccB+$wamaccA+2)." ".($worder+1);
    system "$cmd_WAMBuilder $arg $wam_acc_tp $wam_acc_tp. >& $wam_acc_tp.log";
    system "$cmd_WAMBuilder $arg $wam_acc_fp $wam_acc_fp. >& $wam_acc_fp.log";
    print "StopTF..";
    $arg = "$worder ".($wamstoB+$wamstoA+3)." ".($worder+1);
    system "$cmd_WAMBuilder $arg $wam_sto_tp $wam_sto_tp. >& $wam_sto_tp.log";
    system "$cmd_WAMBuilder $arg $wam_sto_fp $wam_sto_fp. >& $wam_sto_fp.log";
    print ".done\n";

    ## Build logos ##
    print "Launch $cmd_seqlogo (Build logo)....";
    $arg =  "-b -c -n -Y -F PNG";
    print "StartTF..";
    system "$cmd_seqlogo $arg -x StartTP -f $wam_sta_tp > $wam_sta_tp.png";
    system "$cmd_seqlogo $arg -x StartFP -f $wam_sta_fp > $wam_sta_fp.png";
    print "DonTF..";
    system "$cmd_seqlogo $arg -x DonorTP -f $wam_don_tp > $wam_don_tp.png";
    system "$cmd_seqlogo $arg -x DonorFP -f $wam_don_fp > $wam_don_fp.png";
    print "AccTF..";
    system "$cmd_seqlogo $arg -x AcceptorTP -f $wam_acc_tp > $wam_acc_tp.png";
    system "$cmd_seqlogo $arg -x AcceptorFP -f $wam_acc_fp > $wam_acc_fp.png";
    print "StopTF..";
    system "$cmd_seqlogo $arg -x StopTP -f $wam_sto_tp > $wam_sto_tp.png";
    system "$cmd_seqlogo $arg -x StopFP -f $wam_sto_fp > $wam_sto_fp.png";
    print ".done\n";
  }

  PrintInfo();
  FileClose();
}


##################
#     SUB        #
##################
#-----------------------------------------------------------------------------#
sub Usage {
  system("pod2text $0");
  exit();
}

#-----------------------------------------------------------------------------#
#-  Verify the options and open the files corresponding to the approach.     -#
#-----------------------------------------------------------------------------#
sub OptVerif {
  # "genomic" option
  if(defined($genomic)) {
    if(!open(GENOMIC, $genomic)) {
      print STDERR "ERROR: Could not open genomic file list $genomic\n";
      exit(2);
    }
  }
  else { &Usage(); }

  # "cDna" "ssCoord" and "coord" options
  if(defined($cDna) && defined($ssCoord) && !defined($coord)) {
    if(!open(CDNA, $cDna)) {
      print STDERR "ERROR: Could not open cDNA file list $cDna\n";
      exit(2);
    }
    if(!open(SSCOORD, $ssCoord)) {
      print STDERR "ERROR: Could not open start-stop list $ssCoord\n";
      exit(2);
    }
  }
  else {
    if(defined($coord) && !defined($cDna) && !defined($ssCoord)) {
      if(!open(COORD_IN, $coord)) {
	print STDERR "ERROR: Could not open coordinate list $coord\n";
	exit(2);
      }
    }
    else { &Usage(); }
  }

  if($offset < 0) { &Usage(); }

  if(defined($outputDir)) {
    mkdir $outputDir;
    mkdir $outputDir."/Genes";
    $coordFile  = $outputDir."/".$coordFile;
    $exonFile   = $outputDir."/".$exonFile;
    $intFile    = $outputDir."/".$intFile;
    $utr5File   = $outputDir."/".$utr5File;
    $utr3File   = $outputDir."/".$utr3File;
    $matFile    = $outputDir."/".$matFile;
    $matLogFile = $outputDir."/".$matLogFile;
    if(defined($wam)) {
      mkdir $outputDir."/WAM";
      $wam_sta_fp = $outputDir."/WAM/".$wam_sta_fp;
      $wam_don_fp = $outputDir."/WAM/".$wam_don_fp;
      $wam_acc_fp = $outputDir."/WAM/".$wam_acc_fp;
      $wam_sto_fp = $outputDir."/WAM/".$wam_sto_fp;
      $wam_sta_tp = $outputDir."/WAM/".$wam_sta_tp;
      $wam_don_tp = $outputDir."/WAM/".$wam_don_tp;
      $wam_acc_tp = $outputDir."/WAM/".$wam_acc_tp;
      $wam_sto_tp = $outputDir."/WAM/".$wam_sto_tp;
      $verif      = $outputDir."/".$verif;
    }
  }
}

#-----------------------------------------------------------------------------#
#-  Close files and remove the temporary files                               -#
#-----------------------------------------------------------------------------#
sub FileClose {
  # Close input files
  close GENOMIC;
  if(defined($cDna)) { close CDNA; close SSCOORD; }
  else               { close COORD_IN; }

  # Close output files
  close EXON; close INTRON; close UTR5; close UTR3;
  if(defined($wam)) {
    close FPSTART; close FPDON; close FPACC; close FPSTOP;
    close TPSTART; close TPDON; close TPACC; close TPSTOP;
    close VERIF;
  }

  # Remove tmp files
  unlink("s2e_t.txt");
  unlink("s2e_t2.txt");
  unlink("fastaREV");
}

#-----------------------------------------------------------------------------#
#-  First approach.                                                          -#
#-----------------------------------------------------------------------------#
sub GenoCoord {
  while (my $genofile=<GENOMIC>) {
    print"------------------------------------".
      "-------------------------------\n";
    print ($nbSeq+1);
    flush STDOUT;
    print STDERR ".";
    $nbSeq++;
    chomp $genofile;

    ## ID extraction ##
    my $id = IDExtract($genofile);
    print " $id :\n";

    ## Read Coord (.gff) file ##
    print "  - Read coordinates file (.gff)...............................";
    my $gfffile = <COORD_IN>;
    chomp $gfffile;
    my @geneCoord = gffUtils::gffGetCoord("$gfffile");
    print "done\n";

    ## Launch extractseq (EMBOSS) : + offset (DNA) ##
    print "  - Launch extractseq (Create new genomic file,".
      " offset:$offset)...";
    my $begin = ( (abs($geneCoord[1]) - $offset > 0) ?
		  abs($geneCoord[1]) - $offset : 1 );
    my $end   = abs($geneCoord[$#geneCoord-1]) + $offset;
    ExtractSeq ($genofile, $begin, $end,
		" $outputDir"."/Genes/".$id.$newdna_ext);
    print "done\n";

    ## Recomputing coord for the new genomic sequence ##
    my $newFastaFile = "$outputDir/Genes/$id$newdna_ext";
    my $sens  = ( ($geneCoord[1] > 0) ? "" : "-");
    if ($sens eq "-") { $nbSeqRev++; }
    my @newCoord = ();

    print "  - Recompute coordinates :\n";
    print "\t-> Write in coord file...........................";
    my $newcoord = 0;
    my $modif    = ( ($begin < $offset) ?
		     $begin-1 : (abs($geneCoord[1])-$offset-1) );
    print COORD_OUT "$id";
    if ($geneCoord[0] == 0) { push(@newCoord, 0); }
    else                    { push(@newCoord, abs($geneCoord[0])-$modif); }
    for (my $i=1; $i<$#geneCoord; $i++) {
      my $c = abs($geneCoord[$i]) - $modif;
      push(@newCoord, "$sens$c");
      print COORD_OUT " $sens$c";
    }
    if ($geneCoord[-1] == 0) { push(@newCoord, 0); }
    else                     { push(@newCoord, abs($geneCoord[-1])-$modif); }
    print COORD_OUT "\n\n";
    print "done\n";
    print "\t-> Write a gff file..............................";
    gffUtils::gffWrite($id, "$outputDir/Genes/$id$gff_ext",
		       "$newFastaFile", @newCoord);
    print "done\n";

    ## Extract Exons - UTR - Introns ##
    print "  - For build IMM matrices :\n";
    ExtractExon2($id, $newFastaFile, @newCoord);
    ExtractUTRs($id, $newFastaFile, $newCoord[0], $newCoord[1],
		$newCoord[-2], $newCoord[-1]);
    ExtractIntr($id, $newFastaFile, $sens, @newCoord);

    ## Build WAM files ##
    if(defined($wam)) {
      print "  - Build the WAM files....";
      WAM_FP_TP($id, $newFastaFile, @newCoord);
      print "..done\n";
    }
  }
}

#-----------------------------------------------------------------------------#
#-  Second approach.                                                         -#
#-----------------------------------------------------------------------------#
# About each step :
#   For each cDNA seq :
#     - Read Start & Stop
#     - Launch the EMBOSS extractseq on the cDNA file (to extract the CDS)
#     - Launch sim4 against the extracted CDS (with options R=2 A=3 P=1 N=1
#       C=14) and generate an output file.
#     - Read in the sim4 output file the gene coordinates (exons-introns).
# !      A problem was found with very small terminal exon. sim4 was not able
# !      to make the correct alignment.
# !      Because of this problem sim2eugene scans all the exon and if one is
# !      minor than $excluSeq :
# !        1) write in the WARNING file the ID of the seq and the exon length
# !        2) the sequence is excluded
#     - Launch the EMBOSS extractseq on the genomic file (to extract the CDS +
#       a maximum offset on each side) and generate a new genomic file
#     - Add in the coord.txt file the gene coordinates (recomputed for the new
#       genomic sequence) and generate a GFF file (ID.fasta.gff).
#     - Extract Intron - Exon - UTR5/3 and add the sequence in each file :
#          -> Launch the EMBOSS extractseq on the cDNA file for exon and
#	      utr5/3
#          -> Launch the EMBOSS extractseq on the entire genomic file for
#             introns, and if the cDNA align by sim4 is the reverse complement
#             launch the EMBOSS revseq on the extracted seq
#-----------------------------------------------------------------------------#
sub GenoSim4 {
  open(WARNING, ">WARNING")   || die "ERROR: Could not open file WARNING\n";

  while (my $cdnafile=<CDNA>) {
    print"------------------------------------".
      "-------------------------------\n";
    print ($nbSeq+1);
    flush STDOUT;
    print STDERR ".";
    my $exclusion = 0;   # short exon ?  ->  Warning file
    $nbSeq++;
    my $genofile = <GENOMIC>;
    chomp $cdnafile;
    chomp $genofile;

    ## ID extraction ##
    my $id = IDExtract($genofile);
    print " $id :\n";

    ## Read Start Stop file ##
    print "  - Read Start and Stop positions..............................";
    my $ssline = <SSCOORD>;
    chomp $ssline;
    my $start;
    my $stop;
    if ($ssline =~ /([0-9]+)\s+([0-9]+)/) {
      $start = $1;
      $stop  = $2;
    }
    print "done\n";

    ## Launch extractseq (EMBOSS) : start -> stop (cDNA) ##
    print "  - Launch extractseq (CDS from cDNA file).....................";
    ExtractSeq ($cdnafile, $start, $stop, "s2e_mRNAStartStop");
    print "done\n";

    my $simOut  = "$outputDir"."/Genes/"."$id"."$sim4_out_ext";
    my $sens  = "";
    my @coord    = ();
    my @newCoord = ();

    ## Launch sim4 ##
    Sim4 ($simOut, $genofile);

    ## Sim4 output parsing ##
    Sim4Parsing ($simOut, \$sens, \@coord, $id, \$exclusion);

    if( $exclusion == 0 ) {
      ## Launch extractseq (EMBOSS) : + offset (DNA) ##
      print "  - Launch extractseq (Create new genomic file,".
	" offset:$offset)...";
      my $newFastaFile = "$outputDir/Genes/$id$newdna_ext";
      my $begin    = ( ($coord[0] - $offset > 0) ? $coord[0] - $offset : 1 );
      my $end      = $coord[$#coord] + $offset;
      ExtractSeq ($genofile, $begin, $end,
		  " $newFastaFile");
      print "done\n";

      ## Recomputing coord for the new genomic sequence ##
      print "  - Recompute coordinates :\n";
      print "\t-> Write in coord file...........................";
      my $modif    = (($begin<$offset) ? $begin-1 : ($coord[0]-$offset-1));
      my $newcoord = 0;

      print COORD_OUT "$id";
      push(@newCoord, 0);
      foreach $newcoord (@coord) {
	my $c = $newcoord - $modif;
	push(@newCoord, "$sens$c");
	print COORD_OUT " $sens$c";
      }
      push(@newCoord, 0);
      print COORD_OUT "\n\n";
      print "done\n";
      print "\t-> Write a gff file..............................";
      gffUtils::gffWrite($id, "$outputDir/Genes/$id$gff_ext",
			 "$newFastaFile", @newCoord);
      print "done\n";

      ## Extract Exons - UTR - Introns ##
      print "  - For build IMM matrices :\n";
      ExtractExon($id, $cdnafile, $start, $stop);
      ExtractUTRs($id, $cdnafile, $start, $start, $stop, 10000000);
      unshift(@coord, 0);
      push(@coord,    0);
      ExtractIntr($id, $genofile, $sens, @coord);

      ## Build WAM files ##
      if(defined($wam)) {
	print "  - Build the WAM files....";
	WAM_FP_TP($id, $newFastaFile, @newCoord);
	print "..done\n";
      }
    }
    else {
      print "\n  ----> EXCLUSION : contains a very short exon (< $excluSeq)";
      print "    ! WARNING !\n";
    }
  }

  close WARNING;
  if ($nbExclu == 0) { unlink("WARNING"); }
  unlink("s2e_mRNAStartStop");
}

#-----------------------------------------------------------------------------#
#-  ID extraction from genomic file. (1) & (2)                               -#
#-----------------------------------------------------------------------------#
sub IDExtract {
  my $lgenofile = shift;
  if ($lgenofile =~ /\/?(.+\/)*([^\.]+)/) { return $2; }
  else { die "ERROR: Incorrect genomic list (line:$lgenofile)\n" };
}

#-----------------------------------------------------------------------------#
#-  Launch extractseq. (1) & (2)                                             -#
#-----------------------------------------------------------------------------#
sub ExtractSeq {
  my ($seq, $begin, $end, $out) = @_;
  my $regions = "-regions $begin-$end";
  my $trash   = " >& /dev/null";
  system "$cmd_extractseq $seq $regions -outseq $out $trash";
}

#-----------------------------------------------------------------------------#
#-  Launch sim4. (2)                                                         -#
#-----------------------------------------------------------------------------#
sub Sim4 {
  my ($lout, $lgeno) = @_;
  print "  - Launch sim4 (CDS against genomic file).....................";
  system "$cmd_sim4 s2e_mRNAStartStop $lgeno R=2 A=3 P=1 N=1 C=14 > $lout";
  print "done\n";
}

#-----------------------------------------------------------------------------#
#-  Parse the sim4 output for extract the coordinates. (2)                   -#
#-----------------------------------------------------------------------------#
sub Sim4Parsing {
  my ($lsimOut, $refsens, $refcoord, $lid, $refexclu) = @_;
  open (SIM,$lsimOut) || die "ERROR: Could not open file $lsimOut\n";

  while (my $line=<SIM>) {
    if ($line =~ /\(complement\)/) { $$refsens = "-"; $nbSeqRev++; }

    if ($line =~ /\(([0-9]+)\-([0-9]+)\)/) {
      push(@{$refcoord}, $1);
      push(@{$refcoord}, $2);
      my $lenEx = $2-$1+1;
      if( $lenEx < $excluSeq ) {
	print WARNING $lid." :\t".$lenEx."\n";
	$$refexclu = 1;
	$nbExclu++;
      }
      if( $lenEx > $exonMaxLen) { $exonMaxLen = $lenEx; }
      if( $lenEx < $exonMinLen) { $exonMinLen = $lenEx; }
    }
  }
}

#-----------------------------------------------------------------------------#
#-  Extract exon (2).                                                        -#
#-----------------------------------------------------------------------------#
sub ExtractExon {
  my ($lid, $lcdnafile, $lstart, $lstop) = @_;
  print "\t-> Launch extractseq (exon from cDNA file).......";
  print EXON "$lid ";
  ExtractSeq ($lcdnafile, $lstart, $lstop, "s2e_t.txt");
  WriteExtract("s2e_t.txt", *EXON);
  print "done\n";
}

#-----------------------------------------------------------------------------#
#-  Extract exon. (1)                                                        -#
#-----------------------------------------------------------------------------#
sub ExtractExon2 {
  my ($lid, $lfastafile, @lcoord) = @_;

  print "\t-> Launch extractseq (exon from cDNA file)....";
  print EXON "$lid ";
  my $regions = "-regions ";
  my $trash   = " >& /dev/null";
  my $b = "";
  my $e = "";
  for (my $i=1; $i<$#lcoord; $i+=2) {
    $b = abs($lcoord[$i]);
    $e = abs($lcoord[$i+1]);
    $regions .= "$b-$e,";
    my $lenEx = $e - $b + 1;
    if( $lenEx > $exonMaxLen) { $exonMaxLen = $lenEx; }
    if( $lenEx < $exonMinLen) { $exonMinLen = $lenEx; }
  }

  system "$cmd_extractseq $lfastafile $regions -outseq s2e_t.txt $trash";

  if ($lcoord[1] < 0) {
    print "R";
    system "$cmd_revseq s2e_t.txt s2e_t2.txt $trash";
    WriteExtract("s2e_t2.txt", *EXON);
    print ".";
  }
  else {
    WriteExtract("s2e_t.txt", *EXON);
    print "..";
  }
  print ".done\n";
}

#-----------------------------------------------------------------------------#
#-  Extract utr5/3. (1) & (2)                                                -#
#-----------------------------------------------------------------------------#
sub ExtractUTRs {
  my ($lid, $lcdnafile, $l5start, $l5stop, $l3start, $l3stop) = @_;

  ## UTR5 ##
  print "\t-> Launch extractseq (UTR5 from cDNA file)....";
  if (abs($l5start)-1 >= 1) {
    print UTR5 "$lid ";
    if(abs($l5start) == abs($l5stop)) { $l5start = 1; }
    ExtractSeq ($lcdnafile, abs($l5start), abs($l5stop)-1, "s2e_t.txt");
    if ($l5stop < 0) {
      print "R";
      system "$cmd_revseq s2e_t.txt s2e_t2.txt >& /dev/null";
      WriteExtract("s2e_t2.txt", *UTR5);
      print "..";
    }
    else {
      WriteExtract("s2e_t.txt", *UTR5);
      print "...";
    }
    print "done\n";
  }
  else {
    print "...no utr5\n";
    $nbNoUTR5++;
  }
  ## UTR3 ##
  print "\t-> Launch extractseq (UTR3 from cDNA file)....";
  if (abs($l3stop) != 0) {
    ExtractSeq ($lcdnafile, abs($l3start)+1, abs($l3stop), "s2e_t.txt");
    if ($l3start < 0) {
      print "R";
      system "$cmd_revseq s2e_t.txt s2e_t2.txt >& /dev/null";
    }
    else { print "."; }
    open (TMP, "s2e_t.txt") || die "ERROR: Could not open file s2e_t.txt\n";
    my @linesTMP = <TMP>;
    close TMP;
    shift(@linesTMP);
    my $line = "$linesTMP[0]";
    if (length($line) > 2) {
      print UTR3 "$lid ";
      foreach $line (@linesTMP) { chomp $line; print UTR3 $line; }
      print UTR3 "\n";
      print "..done\n";
    }
    else {
      print "..no utr3\n";
      $nbNoUTR3++;
    }
  }
  else {
    print "...no utr3\n";
    $nbNoUTR3++;
  }
}

#-----------------------------------------------------------------------------#
#-  Extract intron. (1) & (2)                                                -#
#-----------------------------------------------------------------------------#
sub ExtractIntr {
  my ($lid, $lgenofile, $lsens, @lcoord) = @_;

  if ($#lcoord != 3) {
    print "\t-> Launch extractseq (intron from genomic file) :\n\t   - ";
    for (my $i=2; $i<$#lcoord-1; $i+=2) {
      print ($i/2);
      print INTRON "$lid.[".($i/2)."] ";
      ExtractSeq ($lgenofile, abs($lcoord[$i])+1,
		  abs($lcoord[$i+1])-1, "s2e_t.txt");
      my $lenInt = abs($lcoord[$i+1])-1 - abs($lcoord[$i]);
      if( $lenInt > $intronMaxLen) { $intronMaxLen = $lenInt; }
      if( $lenInt < $intronMinLen) { $intronMinLen = $lenInt; }
      if ($lsens eq "-") {
	print "R";
	system "$cmd_revseq s2e_t.txt s2e_t2.txt >& /dev/null";
	WriteExtract("s2e_t2.txt", *INTRON);
	print ".";
      }
      else {
	WriteExtract("s2e_t.txt", *INTRON);
	print ".";
      }
    }
    print ".done\n";
  }
  else {
    print "\t-> Launch extractseq (intron from genomic file)..";
    print "no intron\n";
    $nbNoIntron++;
  }
}

#-----------------------------------------------------------------------------#
#-  Write the extracted sequences (exon, intron, utr5/3) in files. (1) & (2) -#
#-----------------------------------------------------------------------------#
sub WriteExtract {
  my ($fileIn, $fileOut) = @_;
  open (TMP, $fileIn) || die "ERROR: Could not open file $fileIn\n";
  my @linesTMP = <TMP>;
  close TMP;
  my $line = "";
  shift(@linesTMP);
  foreach $line (@linesTMP) { chomp $line; print $fileOut $line; }
  print $fileOut "\n";
}

#-----------------------------------------------------------------------------#
#-  Extract and write signals with context for WAMBuilder. (1) & (2)         -#
#-----------------------------------------------------------------------------#
sub WAM_FP_TP {
  my ($lid, $fasta, @coord) = @_;
  my $cmp = 1;
  my $pos;
  my $seq = "";
  my $start = $coord[1];
  my $stop  = $coord[-2];
  my @don = ();
  my @acc = ();
  for (my $i=2; $i<$#coord-1; $i+=2) {
    push(@don, $coord[$i]);
    push(@acc, $coord[$i+1]);
  }

  if($coord[1] < 0) {
    system "$cmd_revseq $fasta fastaREV >& /dev/null";
    $fasta = "fastaREV";
  }

  open(FASTA,$fasta) || die "ERROR: Could not open file $fasta\n";
  my @fasta = <FASTA>;
  close FASTA;
  shift(@fasta);
  foreach my $line (@fasta) { chomp $line; $seq .= $line; }
  my $seqLen = length($seq);
  $seq =~ tr/atcg/ATCG/;

  if($coord[1] < 0) {
    my @coordREV = ();
    foreach my $c (@coord) { unshift(@coordREV, $seqLen - abs($c) + 1); }
    $start = $coordREV[1];
    $stop  = $coordREV[-2];
    @don = ();
    @acc = ();
    for (my $i=2; $i<$#coordREV-1; $i+=2) {
      push(@don, $coordREV[$i]);
      push(@acc, $coordREV[$i+1]);
    }
  }

  print VERIF "$lid\t";

  ## START ##
  my $v = "ATG:NO";
  my $before = $wamstaB + $worder;
  print "Start...";
  while ($seq =~ /(?=([ATCG]{$before}(ATG)[ATCG]{$wamstaA}))/g) {
    $pos = pos($seq) + $before + 1;
    if($pos != $start) {
      print FPSTART ">$lid ".$cmp++." ATG:$pos\n$1\n";
    }
    else {
      print TPSTART ">$lid ".$cmp++." ATG:$pos\n$1\n";
      $v = "ATG:OK";
    }
  }
  print VERIF "$v ";

  ## DONOR ##
  $v = 0;
  if(@don == 0) { push(@don, -1); }
  $cmp  = 1;
  my $j = 0;
  my $donor = $don[$j];
  $before = $wamdonB + $worder;
  print "Donor...";
  while ($seq =~ /(?=([ATCG]{$before}(GT)[ATCG]{$wamdonA}))/g) {
    $pos = pos($seq) + $before;
    if($j < $#don && $pos > $donor) { $donor = $don[++$j]; }
    if($pos != $donor) {
      print FPDON ">$lid ".$cmp++." GT:$pos\n$1\n";
    }
    else {
      print TPDON ">$lid ".$cmp++." GT:$pos\n$1\n";
      $v++;
    }
  }
  if($don[0] != -1) {
    print VERIF " GT$v/".($#don+1);
    ( ($v/($#don+1) == 1) ? print VERIF ":OK " : print VERIF ":NO ");
  }

  ## ACCEPTOR ##
  $v = 0;
  if(@acc == 0) { push(@acc, -1); }
  $cmp = 1;
  $j = 0;
  my $acceptor = $acc[$j];
  $before = $wamaccB + $worder;
  print "Acceptor...";
  while ($seq =~ /(?=([ATCG]{$before}(AG)[ATCG]{$wamaccA}))/g) {
    $pos = pos($seq) + $before + 1 + 2;
    if($j < $#acc && $pos > $acceptor) { $acceptor = $acc[++$j]; }
    if($pos != $acceptor) {
      print FPACC ">$lid ".$cmp++." AG:$pos\n$1\n";
    }
    else {
      print TPACC ">$lid ".$cmp++." AG:$pos\n$1\n";
      $v++;
    }
  }
  if($acc[0] != -1) {
    print VERIF " AG$v/".($#acc+1);
    ( ($v/($#acc+1) == 1) ? print VERIF ":OK " : print VERIF ":NO ");
  }

  ## STOP ##
  $v = "Stop:NO";
  $cmp = 1;
  $before = $wamstoB + $worder;
  print "Stop...";
  my @st  = ("TAA", "TAG", "TGA");
  for(my $i=0; $i<3; $i++) {
    while($seq =~ /(?=([ATCG]{$before}($st[$i])[ATCG]{$wamstoA}))/g) {
      $pos = pos($seq) + $before + 1 + 2;
      if($pos != $stop) {
	print FPSTOP ">$lid ".$cmp++." $st[$i]:$pos\n$1\n";
      }
      else {
	print TPSTOP ">$lid ".$cmp++." $st[$i]:$pos\n$1\n";
	$v = "$st[$i]:OK";
      }
    }
  }
  print VERIF " $v\n";
}

#-----------------------------------------------------------------------------#
#-  Print (stdout) some informations about sequences. (1) & (2)              -#
#-----------------------------------------------------------------------------#
sub PrintInfo {
  print "----------------------------------".
    "---------------------------------\n";
  if ($nbExclu != 0) {
    print"! WARNING $nbExclu / $nbSeq sequence(s) have been excluded".
      " because of very\n!         short exon(s) ( < $excluSeq pb ).\n".
	"!         You MUST read the WARNING file.\n";
    print"---------------------------------".
      "----------------------------------\n";
  }
  print "Informations (for $nbSeq sequence(s)) :\n";
  print "   - Number of :\n";
  print "\t-> Sequences forward : ".($nbSeq-$nbSeqRev)."\n";
  print "\t-> Sequences reverse : $nbSeqRev\n";
  print "\t-> Sequences without utr5   : $nbNoUTR5\n";
  print "\t-> Sequences without utr3   : $nbNoUTR3\n";
  print "\t-> Sequences without intron : $nbNoIntron\n";
  print "   - Length of :\n";
  print "\t-> Exon min.   : $exonMinLen\n";
  print "\t-> Exon max.   : $exonMaxLen\n";
  print "\t-> Intron min. : $intronMinLen\n";
  print "\t-> Intron max. : $intronMaxLen\n";
  system('echo "Finished on "`date`');
  print "-----------------------------------".
    "--------------------------------\n";
  flush STDOUT;
  print STDERR "\n";
}
