//-----------------------------------------------------------------------------
//  This program finds exons/introns and intergenic regions (including UTR)
//  Copyright T. Schiex 1999
//
//  $Id$
// 
// Wed Mar  3 12:06:19 1999 : des starts factices en extremites de seq
// Wed Mar  3 12:28:08 1999 : Indep_eval decale de 1 
// Mon Mar  8 13:27:04 1999 : abaisse la force des START (scale)
// Mon Mar  8 17:27:12 1999 : supprime s/ es. Data commence a 0
// Mon Mar 15 10:12:30 1999 : passage a une version eucaryote
// Mon Mar 15 15:37:00 1999 : patch sigsegv annick - delete Data => (Data-1)
// Thu Mar 18 23:01:35 1999 : usage de macs decale
// Thu Mar 18 23:03:35 1999 : prise en compte de HITS BLAST EST/PROTEINES
// Sun May 23 15:30:26 1999 : modele intronique
// Tue Jun 22 16:50:21 1999 : modele intronique trop fort. 
//                            Utilise en intergenique too
// Tue Jul 6  13:01:23 1999 : modele intronique dirige,
//                            modele intergenique non.
// Mon Jul 19 16:07:54 1999 : introduction de StartB et separation
//                            des parametres Acc et Don
// Mon Jul 19 16:09:54 1999 : Le Start est traite comme une proba cond.
// Mon Oct 18 14:58:36 1999 : macs => backP, traitement des longueurs min

// Fri Mar 16 18:07:10 2001: rel 1.1a: separation du traitement des
//                             EST et correction BestP sous dimensionne

// TODO:
// supprimer Choice
// remettre Frameshifts
// alignement cds cDNA et proteines sur les splice sites

// supprimer l'approximation longueur single (on peut rater un bon
// epissage si un meilleur START l'occulte mais qu'il est trop pres
// pour etre utilisable). Il faudrait faire 3 pistes single + 3 pistes
// non single.

#define  VERSION "1.1d (120701)"


#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <limits.h>
#include <float.h>
#include <time.h>
#include <assert.h>
#include <errno.h>
#ifdef HAVE_GETOPT_H
#include <getopt.h>
#else
#include "getopt.h"
#endif
#include <string.h>
#ifdef HAVE_STRINGS_H
#include <strings.h>
#endif
#include "Const.h"

//#include "dmalloc.h"
//#include "dmallocc.cc"
//#define STAND 1

#define DFT_MATRIX          "default.mat"

#ifndef FILENAME_MAX
#define FILENAME_MAX        1024
#endif



// ------------------ Globals ---------------------
// Les Hits EST

const unsigned char HitForward    = 0x1;
const unsigned char MLeftForward  = 0x2;
const unsigned char GapForward    = 0x4;
const unsigned char MRightForward = 0x8; 

const unsigned char HitReverse    = 0x10;
const unsigned char MLeftReverse  = 0x20;
const unsigned char GapReverse    = 0x40;
const unsigned char MRightReverse = 0x80; 

const unsigned char Hit           = HitForward    | HitReverse;
const unsigned char MLeft         = MLeftForward  | MLeftReverse;
const unsigned char Gap           = GapForward    | GapReverse;
const unsigned char MRight        = MRightForward | MRightReverse;
const unsigned char Margin        = MRight        | MLeft;

const  unsigned char NotAHit       = MLeft | MRight | Gap;

#define Inconsistent(x) (((x) & Hit) && ((x) & Gap))


double Score [9],NScore[9];
int Data_Len;
int normopt,blastopt,estopt,ncopt;
int window, offset,graph,resx,resy;
int    gfrom,gto,golap,glen;

extern char     *optarg;   
extern int      optind;

#include "System.h"
#include "DNASeq.h"
#include "Stop.h"
#include "BStrArray.h"
#include "EuIMMScoreAS.h"
#include "EuStart.h"
#include "BackP.h"
#include "Plot.h"
#include "Hits.h"

//#include "SpliceGS.h"
//#include "SpliceP.h"
//#include "SpliceN.h"
#include "SpliceNP.h"

extern "C"{
#include "../GDIF/gdIF.h"
}
#include "clef.h"


// ------------------ Globals ---------------------

const unsigned char SwitchMask[18] = 
  {SwitchAny,SwitchAny,SwitchAny,SwitchAny,SwitchAny,SwitchAny,
   SwitchAny,SwitchAny,SwitchAny,SwitchAny,SwitchAny,SwitchAny,
   SwitchStart,SwitchStart,SwitchStart,SwitchStop,SwitchStop,SwitchStop};

//--------------------------------------------------
int HitsCompare(const void *A, const void *B)
{
  Hits **UA,**UB;
  
  UA = (Hits **) A;
  UB = (Hits **) B;

  if ((*UA)->NGaps > (*UB)->NGaps) return -1;
  if ((*UA)->NGaps < (*UB)->NGaps) return 1;
  if ((*UA)->Length > (*UB)->Length) return -1;
  if ((*UA)->Length < (*UB)->Length) return 1;
  return 0;
}


int IsPhaseOn(char p, int pp)
{
  if (p == 6) return FALSE;
  if (p == pp) return TRUE;
  return FALSE;
}


// Convertit les phases 0-6 en 1 2 3 -1 -2 -3 0
int PhaseAdapt(char p)
{
  if (p >= 12) return 0;
  else if (p < 3) return (1+p);
  else if (p < 6) return (2-p);
  else if (p < 9) return (p-2);
  else return (5-p);
  
}

void PrintPhase(char p)
{
  switch (p) {
  case 0:
    printf("E1 ");
    return;

  case 1:
    printf("E2 ");
    return;

  case 2:
    printf("E3 ");
    return;

  case 3:
    printf("e1 ");
    return;

  case 4:
    printf("e2 ");
    return;

  case 5:
    printf("e3 ");
    return;

  case 6:
    printf("I1 ");
    return;

  case 7:
    printf("I2 ");
    return;

  case 8:
    printf("I3 ");
    return;

  case 9:
    printf("i1 ");
    return;

  case 10:
    printf("i2 ");
    return;

  case 11:
    printf("i3 ");
    return;

  case 12:
    printf("IG ");
    return;

  case 13:
    printf("U5 ");
    return;

  case 14:
    printf("U3 ");
    return;

  case 15:
    printf("u5 ");
    return;

  case 16:
    printf("u3 ");
    return;
  }
  return;
}

// Convertit les phases 1 2 3 -1 -2 -3 0 en 0-6
char ph06(char p)
{
  if (p == 0) return 6;
  else if (p > 0) return (p-1);
  else return 2-p;   
}

void ESTAnalyzer(FILE *ESTFile,   unsigned char    *ESTMatch,
		 int EstM, REAL** Don, REAL **Acc)
{
  int i,j,k;
  int deb,fin,brin,EstDeb,EstFin,poids,NumEST = 0,Rejected = 0;
  char A[128],B[128];
  char *EstId, *PEstId,*tmp;
  Hits *ThisEST = NULL,*OneEST,*AllEST = NULL;
  Block *ThisBlock = NULL;

  A[0] = B[0]= 0;
  EstId = A;
  PEstId = B;
  
  fprintf(stderr,"Reading cDNA hits...");
  fflush(stderr);

  while (fscanf(ESTFile,"%d %d %d %*s %d %s %d %d\n",
		&deb, &fin,&poids,&brin,EstId,&EstDeb,&EstFin) != EOF)
    {
      if ((strcmp(EstId,PEstId) == 0) && 
	  (EstDeb > ThisBlock->LEnd) && 
	  (deb > ThisBlock->End))
	// si EstId et PEstId sont égaux, alors il y a un EST en cours
	// de meme nom on verifie que c'est bien compatible en terme
	// de position (sur l'est et le genomique)
	  {
	    ThisEST->NGaps++;
	    ThisBlock->AddBlockAfter(deb-1,fin-1,EstDeb,EstFin);
	    ThisBlock = ThisBlock->Next;
	  }
      else {
	NumEST++;
	OneEST = new Hits(EstId,poids,brin,deb-1,fin-1,EstDeb,EstFin);
	ThisBlock = OneEST->Match;
	tmp = PEstId;
	PEstId = EstId;
	EstId = tmp;

	if (AllEST == NULL) {
	  AllEST = OneEST;
	  ThisEST = OneEST;
	}
	else {
	  ThisEST->Next = OneEST;
	  ThisEST = OneEST;
	}
      }
    }
  
  fprintf(stderr,"%d sequences read.\n",NumEST);
  fflush(stderr);

  // on trie les hits sur le nombre de gaps et la
  // longueur. L'idee est d'eliminer les epissages partiels et de
  // favoriser la longueur de toute facon.
  
  Hits **HitTable = new Hits *[NumEST];
  
  for (i = 0, ThisEST = AllEST; i < NumEST; i++, ThisEST = ThisEST->Next)
    HitTable[i] = ThisEST;
  
  qsort((void *)HitTable,NumEST,sizeof(void *),HitsCompare);
  
  for (int index = 0; index < NumEST; index++) {
    
    int Inc;
    REAL WorstSpliceF, WorstSpliceR;
    REAL DonF,AccF,DonR,AccR;
    
    ThisEST = HitTable[index];
    
    Inc = 0;
    WorstSpliceF = WorstSpliceR = 1.0;
    
    
    // Look for each match in the current Hit
    ThisBlock = ThisEST->Match;
    
    while (ThisBlock) {
      // Check for consistency with already read Hits
      
      for (i = ThisBlock->Start+EstM; !Inc && i <= ThisBlock->End-EstM; i++)
	if (Inconsistent(ESTMatch[i] | (HitForward << (ThisEST->Strand*4)))) {
	  fprintf(stderr,"   [%s]: inconsistent hit [%d-%d]\n",
		  ThisEST->Name,ThisBlock->Start+1,ThisBlock->End+1);
	  Inc = 1;
	}
      
      // si on a un gap
      if ((ThisBlock->Prev != NULL) && abs(ThisBlock->Prev->LEnd - ThisBlock->LStart) <= 6) {
	
	for (i = ThisBlock->Prev->End+1+EstM; !Inc && i < ThisBlock->Start-EstM; i++) 
	  if (Inconsistent(ESTMatch[i] | (GapForward << (ThisEST->Strand*4)))) {
	    fprintf(stderr,"   [%s]: inconsistent gap [%d-%d]\n",
		    ThisEST->Name,ThisBlock->Prev->End+2,ThisBlock->Start);
	    Inc = 1;
	  }
	
	// Verification de la presence de site de splicing aux frontieres
	
	DonF = 0.0;
	DonR = 0.0;
	AccF = 0.0;
	AccR = 0.0;
	
	for (j = -EstM; j <= EstM; j++) {
	  k = Min(Data_Len,Max(0,ThisBlock->Prev->End+j+1));
	  DonF = Max(DonF, Don[0][k]);
	  AccR = Max(AccR,Acc[1][k]);
	  
	  k = Min(Data_Len,Max(0,ThisBlock->Start+j));
	  DonR = Max(DonR,Don[1][k]);
	  AccF = Max(AccF, Acc[0][k]);
	}
	WorstSpliceF = Min(WorstSpliceF,DonF);
	WorstSpliceF = Min(WorstSpliceF,AccF);
	WorstSpliceR = Min(WorstSpliceR,DonR);
	WorstSpliceR = Min(WorstSpliceR,AccR);
      }
      
      ThisBlock = ThisBlock->Next;
    }
    
    // Tous les blocs ont ete traites
    
    if ((WorstSpliceF == 0.0) && (WorstSpliceR == 0.0)) {
      fprintf(stderr, "   [%s]: no matching splice site\n",ThisEST->Name);
      Inc =2;
    }
    
    // Si une incoherence est detectee, on va jeter la sequence
    
    if (Inc) {
      Rejected++;
      
      ThisBlock = ThisEST->Match;
      while (ThisBlock) {
	for (i = ThisBlock->Start; i<= ThisBlock->End; i++) {
	  PlotBarI(i, 4,0.9,1,7);
	  PlotBarI(i,-4,0.9,1,7);
	}
	
	if (graph && (ThisBlock->Prev != NULL) && 
	    abs(ThisBlock->Prev->LEnd-ThisBlock->LStart) <= 6) {
	  PlotLine(ThisBlock->Prev->End,ThisBlock->Start,4,4,0.9,0.9,7);
	  PlotLine(ThisBlock->Prev->End,ThisBlock->Start,-4,-4,0.9,0.9,7);
	}
	ThisBlock = ThisBlock->Next;
      }
    }
    // sinon on l'exploite
    else {
      int LBoundary;
      int RBoundary;
      
      ThisBlock = ThisEST->Match;
      
      while (ThisBlock) {
	
	// Aligners tend to extend beyond the true hit on
	// extremities: we remove EstM on frontiers
	
	LBoundary = ((ThisBlock->Prev == NULL) ? 
		     Min(Data_Len, ThisBlock->Start+EstM) :
		     ThisBlock->Start);
	
	RBoundary = ((ThisBlock->Next == NULL) ? 
		     Max(0,ThisBlock->End-EstM) :
		     ThisBlock->End);
	
	for (i = LBoundary; i <= RBoundary; i++)
	  ESTMatch[i] |= (HitForward << (ThisEST->Strand*4));
	
	if ((ThisBlock->Prev != NULL) && 
	    abs(ThisBlock->Prev->LEnd-ThisBlock->LStart) <= 6) {
	  
	  for (i = Max(0,ThisBlock->Prev->End+1-EstM); 
	       i < Min(Data_Len, ThisBlock->Prev->End+1+EstM); i++) 
	    ESTMatch[i] |= (MLeftForward << (ThisEST->Strand*4));
	  
	  for (i = ThisBlock->Prev->End+1; i < ThisBlock->Start-1; i++)
	    ESTMatch[i] |= (GapForward << (ThisEST->Strand*4));
	  
	  for (i =  Max(0,ThisBlock->Start-1-EstM);
	       i <  Min(Data_Len, ThisBlock->Start-1+EstM); i++)
	    ESTMatch[i] |= (MRightForward << (ThisEST->Strand*4));
	  
	  if (graph) {
	    PlotLine(ThisBlock->Prev->End,ThisBlock->Start,4,4,0.6,0.6,2);
	    PlotLine(ThisBlock->Prev->End,ThisBlock->Start,-4,-4,0.6,0.6,2);
	  }
	  
	}
	ThisBlock = ThisBlock->Next;
      }
    }
  }
  
  for (i = 0; i < Data_Len; i++) 
    if (graph && (ESTMatch[i] & Hit)) {
      PlotBarI(i, 4,0.6,1,2);
      PlotBarI(i, -4,0.6,1,2);
    }
  if (Rejected)
    fprintf(stderr,"A total of %d/%d sequences rejected\n",Rejected,NumEST);
  delete AllEST;
  delete HitTable;
}

// -------------------------------------------------------------------------
//            MAIN
// -------------------------------------------------------------------------
int main  (int argc, char * argv [])
{
  int              i, j, k, carg, errflag;
  FILE             *fp;
  REAL             *BaseScore[11],*Start[2] ;
  REAL             *Don[2],*Acc[2];
  unsigned char    *Stop[2];
  unsigned char    *ESTMatch;
  BString_Array    *IMMatrix[6];
  char             clef[20];
  DNASeq           *TheSeq;
  char             printopt;
  char             matname[FILENAME_MAX+1], tempname[FILENAME_MAX+1], fstname[FILENAME_MAX+1];
  char             grname[FILENAME_MAX+1];
  int              EstM;
  double           FsP,StartP,StartB,StopP,TransStartP,TransStopP;
  double           AccP,AccB,DonP,DonB,BlastS[3],EstP;
  unsigned char *ForcedIG;

  char *EugDir;

  // prior on the initial state, selon Sato et al 1999 / Terryn et
  // al. 1999
  REAL ExonPrior = 0.33,IntronPrior = 0.17;
  REAL InterPrior = 0.4,FivePrimePrior = 0.03,ThreePrimePrior = 0.07;

  const REAL DontCrossStop = NINFINITY;
  const REAL IGPenalty = NINFINITY; 

  // Les longueurs. 

  int MinFivePrime, MinThreePrime;
  int MinEx, MinIn, MinSg,MinFlow,MinConv,MinDiv;
  int MinLength[18];

  // process args 
  // default values

  glen = -1;
  golap = -1;
  gfrom = -1;
  gto = -1;
  resx = 900;                 // x res for graph. output
  resy = 400;                 // y res for graph. output
  graph = FALSE;              // don't produce a graphical output
  estopt = FALSE;             // don't try to read a EST file
  blastopt = -1;              // don't try to read a blast file
  ncopt = FALSE;                   //don't try to read a non coding input
  normopt = 1;                // normalize across frames
  window  = 97;               // window length
  printopt = 'l';             // short print format 
  offset = 0;                 // no offset

  (void) strcpy(matname, DFT_MATRIX); // default matrix file
  *fstname = '\000';                  // no default input    
  errflag = 0;


  EugDir = getenv("EUGENEDIR");

  strcpy(tempname,argv[0]);
  strcat(tempname,".par");
  fp = FileOpen(EugDir,BaseName(tempname),"r");

  fprintf(stderr,"EuGene rel. %s\n",VERSION);
  fprintf(stderr,"Loading parameters file...");
  fflush(stderr);
  
  if (fscanf(fp,
"%s %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %d %d %d %d %d %d %d %d %d",
	     clef,&FsP,&StartP,&StartB,&StopP,&TransStartP,&TransStopP,&AccP,&AccB,
	     &DonP,&DonB,&BlastS[0],&BlastS[1],&BlastS[2],&EstP,&EstM,
	     &MinEx,&MinIn,&MinSg,&MinFlow,&MinConv,&MinDiv,
	     &MinFivePrime,&MinThreePrime) != 24)
    {
      fprintf(stderr, "Incorrect parameter file EuGene.par\n");
      exit(2);
    }
  fprintf(stderr,"done\n");
  fclose(fp);

 // remplir le tableau des longueurs min de chaque etat (+ 6 pour les Single)

  for (i = 0; i < 6; i++) {
    MinLength[i] = MinEx;
    MinLength[i+6] = MinIn;
    MinLength[i+12] = MinSg;
  }

  ReadKey(clef,"EUGENEAT");

  // any Frameshift prob below -1000.0 means "not possible"
  if (FsP <= -1000.0) FsP = NINFINITY;

  while ((carg = getopt(argc, argv, "drshm:w:f:n:o:p:x:y:u:v:g::b::l:")) != -1) {
    
    switch (carg) {
            
    case 'n':            /* -n normalize across frames      */
      if (! GetIArg(optarg, &normopt, normopt))
	errflag++;
      break;
      
    case 'h':           /* help                             */
      errflag++;
      break;

    case 's':           /* Single gene mode. Start/End on UTR or IG only   */
      ExonPrior = 0.0;
      IntronPrior = 0.0;
      InterPrior = 1.0;
      FivePrimePrior = 1.0;
      ThreePrimePrior = 1.0;
      break;

    case 'r':    /* RepeatMasker input */
      ncopt = TRUE;
      break;
      
    case 'c':            /* -c couverture      */
      if (! GetIArg(optarg, &golap, golap))
        errflag++;
      break;
      
    case 'l':            /* -l imglen      */
      if (! GetIArg(optarg, &glen, glen))
        errflag++;
      break;
      
    case 'u':            /* -u From      */
      if (! GetIArg(optarg, &gfrom, gfrom))
        errflag++;
      break;
      
    case 'v':            /* -v To      */
      if (! GetIArg(optarg, &gto, gto))
        errflag++;
      break;
      
    case 'x':            /* -x resx      */
      if (! GetIArg(optarg, &resx, resx))
        errflag++;
      break;
      
    case 'y':            /* -y resy      */
      if (! GetIArg(optarg, &resy, resy))
        errflag++;
      break;
      
    case 'g':           /* -g "Graphic File"                */
      if (optarg) strcpy(grname,optarg);
      else grname[0] = 0;
      graph = TRUE;
      break;
    
    case 'p':           /* print opt: short/long/detailed   */
      printopt = optarg[0];
      if ((printopt == 'h') &&  graph == FALSE) {
	graph = TRUE; // HTML output means graphical output 
	grname[0] = 0;
      }
      if ((printopt != 's') && (printopt != 'l') && 
	  (printopt != 'd')&& (printopt != 'h'))
	errflag++;
      break;
      
    case 'm':           /* -m matrix                        */
      (void) strcpy(matname, optarg);
      break;

    case 'o':           /* -o offset                        */
      if (! GetIArg(optarg, &offset, offset))
	errflag++;
      break;
      
    case 'w':           /* -w window                        */
      if (! GetIArg(optarg, &window, window))
	errflag++;
      if (window == 0) window = 1;
      break;

    case 'b':           /* -b  use blastx result            */
      if (optarg)
	errflag += !GetIArg(optarg, &blastopt, blastopt);
      else
	blastopt = 7;

      if (blastopt > 7) errflag++;
      break; 
      
    case 'd':           /* -d  use cDNA blastn results      */
      estopt = TRUE;
      break; 

    case 'f':           /* -f frameshift loglike            */
      if (! GetDArg(optarg, &FsP, FsP))
	errflag++;
      break;
      
    case '?':           /* bad option                       */
      errflag++;
    }
  }


  // may remain arguments -> fasta filenames

  if ((argc - optind) == 0)
    errflag++;
    
  // check usage
  
  if (errflag) {
    fprintf(stderr, "usage: EuGene [-h] [-m matrix] [-n 0|1|2] [-s] [-p h|s|l|d] [-w window]\n");
    fprintf(stderr, "              [-b {level}] [-d] [-o offset] [-g {groutput}] [-u start]\n");
    fprintf(stderr, "              [-v end] [-l len] [-c olap] [-x xres] [-y yres] FASTA files\n");
    exit(1);
  }
  
  // open matrix
    
  if (! (fp = FileOpen(EugDir,matname,  "rb"))) {
    fprintf(stderr, "cannot open matrix file %s\n",  matname);
    exit(2);
  }

  fprintf(stderr,"Loading Markov model...");
  fflush(stderr);

  // On essaie d'abord de charger les 5 modeles fondamentaux (3 exons/introns/interG)
  for  (i = 0;  i < 5;  i ++) {
    IMMatrix[i] = new BString_Array(MODEL_LEN, ALPHABET_SIZE);
    if (IMMatrix[i]->Read(fp)) {
      fprintf(stderr,"Model %d unreadable in %s. Aborting.\n",i+1,matname);
      exit(1);
    }
    
    fprintf(stderr,"%d ",i+1);
    fflush(stderr);
  }
  // On essaie ensuite de lire un 6eme modele UTR3'. Si cela echoue,
  // le modele interG est utilise pour l'UTR 3'

    IMMatrix[5] = new BString_Array(MODEL_LEN, ALPHABET_SIZE);
    if (IMMatrix[i]->Read(fp)) {
      fprintf(stderr,"- No 3'UTR model found, using intergenic model. ");
      delete IMMatrix[5];
      IMMatrix[5] = IMMatrix[4];
    } else
      fprintf(stderr,"6 ");

  fprintf(stderr,"done\n");
  fclose(fp);

  int sequence;
  for (sequence = optind; sequence < argc ; sequence++) {
      
    (void) strcpy(fstname, argv[sequence]);
    
    // read fasta file
    
  fp = (*fstname ? FileOpen (NULL,fstname, "r") : stdin);

  if (fp == NULL) {
    fprintf(stderr, "cannot open fasta file %s\n",  fstname);
    exit(3);
  }

  fprintf(stderr,"Loading sequence...");
  fflush(stderr);
  
  TheSeq = new DNASeq(fp);

  if (fp != stdin) fclose(fp);
  
  Data_Len = TheSeq->SeqLen;

  fprintf(stderr,"%s, %d bases read\n",TheSeq->Name, Data_Len);

  fprintf (stderr,"GC Proportion = %.1f%%\n", 
	   (TheSeq->Markov0[BitG]+TheSeq->Markov0[BitC])*100.0);

  for (i = 0;  i < 11;  i ++)
    BaseScore[i] = new REAL[Data_Len+1];
 
  Stop[0] = new unsigned char[Data_Len+1];
  Stop[1] = new unsigned char[Data_Len+1];

  ESTMatch = new unsigned char[Data_Len+1];

  for (i = 0; i<= Data_Len; i++) ESTMatch[i] = 0;

  Start[0] = new REAL[Data_Len+1];
  Start[1] = new REAL[Data_Len+1];

  Acc[0] = new REAL[Data_Len+1];
  Acc[1] = new REAL[Data_Len+1];

  Don[0] = new REAL[Data_Len+1];
  Don[1] = new REAL[Data_Len+1];
  
  Find_Stop_Codons(TheSeq, Stop);

  fprintf(stderr, "Reading start files...");
  fflush(stderr);
  
  strcpy(tempname,fstname);
  strcat(tempname,".starts");
  Read_Start(tempname,Data_Len, exp(-StartP), StartB, Start[0],0);
  fprintf(stderr,"forward,");
  fflush(stderr);

  strcpy(tempname,fstname);
  strcat(tempname,".startsR");
  Read_Start(tempname,Data_Len, exp(-StartP), StartB, Start[1],1);
  fprintf(stderr," reverse done\n");

  fprintf(stderr, "Reading splice site files...");  
  fflush(stderr);

  Read_Splice(fstname,1,Data_Len, Acc[0], Don[0],AccP,AccB,DonP,DonB);
  fprintf(stderr,"forward,");
  fflush(stderr);

  Read_Splice(fstname,-1,Data_Len, Acc[1], Don[1],AccP,AccB,DonP,DonB);
  fprintf(stderr," reverse done\n");

  if (graph) {
    if (grname[0] == 0) {
      strcpy(grname,fstname);
      *rindex(grname,'.') = 0; // on enleve l'extension (.fasta typ.)
    }

    if (gfrom != -1) sprintf(grname+strlen(grname),"%d",gfrom);
    if (gto != -1) sprintf(grname+strlen(grname),"-%d",gto);
    
    if ((gfrom <= 0)|| (gfrom >= Data_Len)) gfrom = 1;
    if ((gto <= 0)  || (gto <= gfrom) || (gto > Data_Len)) gto = Data_Len;
    
    gfrom--;
    gto--;
    
    if (glen < 0) glen = ((gto-gfrom+1 <= 6000) ? gto-gfrom+1 : 6000);
    
    InitPNG(resx,resy,offset,gfrom,gto,golap,glen,grname);
  }
    
  // Precompute all probabilities
  
  fprintf(stderr, "Computing coding probabilities...");  
  fflush(stderr);
  
  Fill_Score(TheSeq,IMMatrix,BaseScore);

  for (j = 0;  j < 11;  j ++)
    BaseScore[j][Data_Len] = 1.0;

  fprintf(stderr,"done\n");

  OpenDoor();

  if (!PorteOuverte && Data_Len > 6000) 
    {
      printf("Valid license key not found ! The software is therefore limited to 6kb long\n");
      printf("sequences. Please contact Thomas.Schiex@toulouse.inra.fr and ask for a FREE\n");
      printf("license key.\n");
      
      exit(2);
    }
  
  /* Exploiting spliced alignements against EST and complete cDNA */

  if (estopt) {
    FILE *fEST;
    strcpy(tempname,fstname);
    strcat(tempname,".est");
    fEST = FileOpen(NULL,tempname, "r");
    ESTAnalyzer(fEST,ESTMatch,EstM,Don,Acc);
    fclose(fEST);
  }


  if (graph) PlotScore(Data_Len,window,normopt,BaseScore);

  // Analysing NC input
  if (ncopt) {
    FILE* ncfile;
    int deb,fin;

    fprintf(stderr,"Reading Intergenic regions... ");
    fflush(stderr);

    strcpy(tempname,fstname);
    strcat(tempname,".ig");
    ncfile = FileOpen(NULL,tempname, "r");

    ForcedIG = new unsigned char[Data_Len+1];
    for (j = 0; j <= Data_Len; j++) ForcedIG[j] = FALSE;

    while (fscanf(ncfile,"%d %d\n", &deb, &fin) != EOF)  {
      deb = Max(1,deb)-1;
      fin = Min(Data_Len,fin)-1;
      for (i = deb; i <= fin; i++) ForcedIG[i] = TRUE;
    }
  }

  // Blastx against protein databases 1st simple version, we simply
  // enhance coding probability according to phase Dangerous (NR
  // contains translation of frameshifted sequences) another
  // possibility would be too forbid introns/intergenic states
  // 3 levels of confidence may be used.
 

  const int LevelColor[3] = {6,7,8}; 
  if (blastopt > 0)
    {
      FILE *fblast;
       int deb,fin,phase,score,Pfin,ProtDeb,ProtFin,PProtFin,PPhase,level;
       char A[128],B[128];
       char *ProtId, *PProtId,*tmp;
       REAL BScore;
       REAL *BlastScore[6];
       unsigned char *InterBlast;

       fprintf(stderr,"Reading Blastx data, level... ");
       fflush(stderr);
       
       for( level = blastopt; level >= 0; level--)
	 {
	   if (blastopt&(1<<level)){
	     strcpy(tempname,fstname);
	     strcat(tempname,".blast");
	     i = strlen(tempname);
	     tempname[i] = '0'+level;
	     tempname[i+1] = 0;
	     
	     fblast = FileOpen(NULL,tempname, "r");
	     
	     if (fblast == NULL) continue;
	     
	     for (i = 0; i < 6; i++)  {
	       BlastScore[i] = new REAL[Data_Len];
	       for (j=0; j<Data_Len; j++) BlastScore[i][j] = 0.0;
	     }
	     
	     InterBlast  = new unsigned char[Data_Len];
	     for (j = 0; j < Data_Len; j++) InterBlast[j] = 0;
	     
	     A[0] = B[0]= 0;
	     ProtId = A;
	     PProtId = B;
	     PProtFin = -10000;
	     
	     while (fscanf(fblast,"%d %d %d %*s %d %s %d %d\n", 
			   &deb, &fin, &score, &phase, ProtId,&ProtDeb,&ProtFin) != EOF) {
	       if (phase < 0) {
		 j = deb;
		 deb = fin;
		 fin = j;
		 j = ProtDeb;
		 ProtDeb = ProtFin;
		 ProtFin = j;
	       }
	       
	       if ((strcmp(ProtId,PProtId) == 0) && (abs(PProtFin-ProtDeb) <= 8)) {
		 // blastx tends to enlarge local hits, we do not reduce gap width (no margin).
		 for (i = Pfin; i < deb; i++) InterBlast[i] = level+1;
		 if (graph) PlotLine(Pfin,deb,PPhase,phase,0.6+(level/8.0),0.6+(level/8.0),LevelColor[level]);
	       }
	       
	       Pfin = fin;
	       tmp = PProtId;
	       PProtId = ProtId;
	       ProtId = tmp;
	       PProtFin = ProtFin;
	       PPhase = phase;
	       
	       phase = ph06(phase);
	       BScore = ((REAL)score)/((REAL)abs(fin-deb));
	       
	       for (i = deb-1; i < fin; i++)  {
		 if (BScore > BlastScore[phase][i]) 
		   BlastScore[phase][i] = BScore;
	       }
	     }

	     fprintf(stderr,"%d ",level);
	     fflush(stderr);
	     
	     for (i = 0; i < 6; i++)  {
	       for (j = 0; j < Data_Len; j++)  {
		 BaseScore[i][j] = 
		   ((BaseScore[i][j]*100.0)+(BlastScore[i][j]*BlastS[level]))/
		   (100.0+(BlastScore[i][j]*BlastS[level]));
		 if (graph && (BlastScore[i][j] > 0.0)) 
		   PlotBarI(j,PhaseAdapt(i),0.6+(level/8.0),1,LevelColor[level]);
	       }
	       delete BlastScore[i];
	     }
	     
	   /*       
	   // Parametre implicite !!!!!
	   int tmp = 0;
	   for (j = 0; j < Data_Len; j++) { 
	     if (InterBlast[j]) {
	       BaseScore[8][j] = (BaseScore[8][j]*100.0)/(100.0);
	       if (!tmp) {
		 printf("%d-",j);
		 tmp = 1;
	       }
	     }
	     else if (tmp) {
	       tmp = 0;
	       printf("%d\n",j-1);
	     }
	   }
	   */
	     delete InterBlast;
	     
	     fclose(fblast);       
	   }
	 }
       fprintf(stderr," done\n");
    }

  // Data allocation for the shortest path with length constraints algorithm
  //
  // 0  codant 1
  // 1  codant 2
  // 2  codant 3
  // 3  codant -1
  // 4  codant -2
  // 5  codant -3
  // 6  Intron frame 0
  // 7  Intron frame 1
  // 8  Intron frame 2
  // 9  Intron frame 0
  // 10 Intron frame 1
  // 11 Intron frame 2  
  // 12 Intrgenique
  // 13 UTR 5' direct
  // 14 UTR 3' direct
  // 15 UTR 5' reverse
  // 16 UTR 3' reverse
  
  char *Choice;
  BackPoint *LBP[18];
  REAL BestU;
  signed char best;
  unsigned char Switch;

  Choice =  new char[Data_Len+2];
  for (i = 0; i < Data_Len +2; i++) Choice[i] = 0;
  
  for (i = 0; i < 18; i++) {
    LBP[i] = new BackPoint(i, -1,0.0);
    LBP[i]->Next = LBP[i];
    LBP[i]->Prev = LBP[i];
  }
  
  const int ExonF1 = 0;
  const int ExonF2 = 1;
  const int ExonF3 = 2;
  const int ExonR1 = 3;
  const int ExonR2 = 4;
  const int ExonR3 = 5;

  const int IntronF1 = 6;
  const int IntronF2 = 7;
  const int IntronF3 = 8;
  const int IntronR1 = 9;
  const int IntronR2 = 10;
  const int IntronR3 = 11;

  const int InterGen5 = 12;
  const int InterGen3 = 17;

  const int UTR5F = 13;
  const int UTR5R = 15;
  const int UTR3F = 14;
  const int UTR3R = 16;

  
  // Codant
  LBP[ExonF1]->Update(log(ExonPrior/6.0)/2.0);
  LBP[ExonF2]->Update(log(ExonPrior/6.0)/2.0);
  LBP[ExonF3]->Update(log(ExonPrior/6.0)/2.0);
  LBP[ExonR1]->Update(log(ExonPrior/6.0)/2.0);
  LBP[ExonR2]->Update(log(ExonPrior/6.0)/2.0);
  LBP[ExonR3]->Update(log(ExonPrior/6.0)/2.0);
  
  // Intron
  LBP[IntronF1]->Update(log(IntronPrior/6.0)/2.0);
  LBP[IntronF2]->Update(log(IntronPrior/6.0)/2.0);
  LBP[IntronF3]->Update(log(IntronPrior/6.0)/2.0);
  LBP[IntronR1]->Update(log(IntronPrior/6.0)/2.0);
  LBP[IntronR2]->Update(log(IntronPrior/6.0)/2.0);
  LBP[IntronR3]->Update(log(IntronPrior/6.0)/2.0);
  
  // Intergenique 
  LBP[InterGen5]->Update(log(InterPrior)/2.0); // ameliorations sur 5' reverse
  LBP[InterGen3]->Update(log(InterPrior)/2.0); // ameliorations sur 3' direct

  // UTR 5' et 3'
  LBP[UTR5F]->Update(log(FivePrimePrior/2.0)/2.0);
  LBP[UTR3F]->Update(log(ThreePrimePrior/2.0)/2.0);  
  LBP[UTR5R]->Update(log(FivePrimePrior/2.0)/2.0);
  LBP[UTR3R]->Update(log(ThreePrimePrior/2.0)/2.0);

  // Les PrevBP sont des pointeurs sur les "opening edges"
  // Les PBest correspondent au cout du chemin correspondant

  REAL  maxi, PBest[26];
  BackPoint *PrevBP[26];
  int source;

  for (i = 0; i <= Data_Len; i++) {
    
    // Calcul des meilleures opening edges

    PrevBP[0] = LBP[0]->BestUsable(i,SwitchMask[0],MinLength[0],&PBest[0]);
    PrevBP[1] = LBP[1]->BestUsable(i,SwitchMask[1],MinLength[1],&PBest[1]);
    PrevBP[2] = LBP[2]->BestUsable(i,SwitchMask[2],MinLength[2],&PBest[2]);

    for (k = 3 ; k < 18; k++) 
      PrevBP[k] = LBP[k%12]->BestUsable(i,SwitchMask[k],MinLength[k],&PBest[k]);
    
    // intergenic: longueur min depend du sens
    // -> -> ou <- <-   MinFlow
    // -> <-            MinConv
    // <- ->            MinDiv
    PrevBP[23] = LBP[InterGen5]->StrictBestUsable(i,MinDiv,&PBest[23]);
    PrevBP[24] = LBP[InterGen5]->StrictBestUsable(i,MinFlow,&PBest[24]);

    PrevBP[18] = LBP[InterGen3]->StrictBestUsable(i,MinFlow,&PBest[18]);
    PrevBP[25] = LBP[InterGen3]->StrictBestUsable(i,MinConv,&PBest[25]);

    // UTR 5' et 3' direct
    PrevBP[19] = LBP[UTR5F]->StrictBestUsable(i,MinFivePrime,&PBest[19]);
    PrevBP[20] = LBP[UTR3F]->StrictBestUsable(i,MinThreePrime,&PBest[20]);
    // UTR 5' et 3' reverse
    PrevBP[21] = LBP[UTR5R]->StrictBestUsable(i,MinFivePrime,&PBest[21]);
    PrevBP[22] = LBP[UTR3R]->StrictBestUsable(i,MinThreePrime,&PBest[22]);

    // ----------------------------------------------------------------
    // ------------------ Exons en forward ----------------------------
    // ----------------------------------------------------------------
    for (k = 0; k < 3; k++) {
      
      maxi = NINFINITY;     

      // S'il y  a un STOP en phase on ne peut continuer
      if ((i % 3 == k) && Stop[0][i]) 
	LBP[k]->Update(DontCrossStop); 
      else // on ne prend pas le donneur
	LBP[k]->Update(log(1.0-Don[0][i]));
      
      LBP[k]->BestUsable(i,SwitchAny,0,&BestU);
      // Un test tordu pour casser le cou aux NaN
      if (isnan(maxi) || (BestU > maxi)) {
	maxi = BestU;
	best = k;
      }
      
      // On commence a coder (Start)
      // Ca vient d'une UTR 5' forward

      if ((i % 3 == k) && Start[0][i] != 0.0) {
	BestU = PBest[19]+log(Start[0][i]);
	// Un test tordu pour casser le cou aux NaN
	if (isnan(maxi) || (BestU > maxi)) {
	  maxi = BestU;
	  best = 19;
	  source = 13;
	  Switch = SwitchStart;
	}
      }
      
      // On recommence a coder (Accepteur)
      // Ca vient d'un intron

      BestU = PBest[6+((i-k+3) % 3)]+log(Acc[0][i]);
      // Un test tordu pour casser le cou aux NaN
      if (isnan(maxi) || (BestU > maxi)) {
	maxi = BestU;
	best = 6+((i-k+3) % 3);
	Switch = SwitchAcc;
      }
      
      if (best != k) 
	LBP[k]->InsertNew(((best >= 18) ? source : best),Switch,i,maxi,PrevBP[best]);
      
      LBP[k]->Update(log(BaseScore[k][i]));
      LBP[k]->Update(((ESTMatch[i] & Gap) != 0)*EstP);

      if ((ForcedIG != NULL) && ForcedIG[i])
	LBP[k]->Update(IGPenalty);
    }
    // ----------------------------------------------------------------
    // ------------------------- Exons en reverse ---------------------
    // ----------------------------------------------------------------
    for (k = 3; k<6; k++) {

      maxi = NINFINITY;

      // On continue sauf si l'on rencontre un autre STOP
      if (((Data_Len-i) % 3 == k-3) && Stop[1][i]) 
	LBP[k]->Update(DontCrossStop);
      else // sinon on ne prend pas le donneur 
	LBP[k]->Update(log(1.0-Don[1][i]));

      
      LBP[k]->BestUsable(i,SwitchAny,0,&BestU);
      // Un test tordu pour casser le cou aux NaN
      if (isnan(maxi) || (BestU > maxi)) {
	maxi = BestU;
	best = k;
      }

      // On commence a coder (Stop)
      // Ca vient d'une UTR 3' reverse

      if (((Data_Len-i) % 3 == k-3) && Stop[1][i]) {
	BestU = PBest[22]-StopP;
	// Un test tordu pour casser le cou aux NaN
	if (isnan(maxi) || (BestU > maxi)) {
	  maxi = BestU;
	  best = 22;
	  source = 16;
	  Switch = SwitchStop;
	}
      }

      // - on recommence a coder (Donneur)
      // Ca vient d'un intron

      BestU = PBest[9+((Data_Len-i-k) % 3)]+log(Don[1][i]);
      // Un test tordu pour casser le cou aux NaN
      if (isnan(maxi) || (BestU > maxi)) {
	maxi = BestU;
	best = 9+((Data_Len-i-k) % 3);
	Switch =  SwitchDon;
      }
            
      if (best != k) 
	LBP[k]->InsertNew(((best >= 19) ? source : best),Switch,i,maxi,PrevBP[best]);
      
      LBP[k]->Update(log(BaseScore[k][i])) ;
      LBP[k]->Update(((ESTMatch[i] & Gap ) != 0)*EstP);

      if ((ForcedIG != NULL) && ForcedIG[i])
	LBP[k]->Update(IGPenalty);
    }
    // ----------------------------------------------------------------
    // ------------------------ Intergenique --------------------------
    // ----------------------------------------------------------------
    // Ca peut venir d'une fin de 3' direct ou de 5' reverse


    // ----------------------------------------------------------------
    // D'abord les 5' reverse => LBP[InterGen5]
    // ----------------------------------------------------------------
    
    maxi = NINFINITY;
    
    // On reste intergenique
    LBP[InterGen5]->BestUsable(i,SwitchAny,0,&BestU);
    // Un test tordu pour casser le cou aux NaN
    if (isnan(maxi) || (BestU > maxi)) {
      maxi = BestU;
      best = -1;
    }

    // From 5' reverse
    BestU = PBest[21]-TransStartP;
    // Un test tordu pour casser le cou aux NaN
    if (isnan(maxi) || (BestU > maxi)) {
      maxi = BestU;
      best = 21;
      source  = 15;
      Switch = SwitchTransStart;
    }
    
    if (best != -1) 
      LBP[InterGen5]->InsertNew(source,Switch,i,maxi,PrevBP[best]);
    
    LBP[InterGen5]->Update(log(BaseScore[8][i]));
    LBP[InterGen5]->Update(((ESTMatch[i] & (Gap|Hit)) != 0)*EstP);

    // ----------------------------------------------------------------
    // Puis les 3' forward => LBP[InterGen3]
    // ----------------------------------------------------------------
    
    maxi = NINFINITY;
    
    // On reste intergenique
    LBP[InterGen3]->BestUsable(i,SwitchAny,0,&BestU);
    // Un test tordu pour casser le cou aux NaN
    if (isnan(maxi) || (BestU > maxi)) {
      maxi = BestU;
      best = -1;
    }
    
    // From 3' direct
    BestU = PBest[20]-TransStopP;
    // Un test tordu pour casser le cou aux NaN
    if (isnan(maxi) || (BestU > maxi)) {
      maxi = BestU;
      best = 20;
      source  = 14;
      Switch = SwitchTransStop;
    }
   
    if (best != -1) 
      LBP[InterGen3]->InsertNew(source,Switch,i,maxi,PrevBP[best]);
    
    LBP[InterGen3]->Update(log(BaseScore[8][i]));
    LBP[InterGen3]->Update(((ESTMatch[i] & (Gap|Hit)) != 0)*EstP);

    // ----------------------------------------------------------------
    // ---------------------- UTR 5' direct ---------------------------
    // ----------------------------------------------------------------
    maxi = NINFINITY;
    
    // On reste 5' direct. On ne prend pas le Start eventuel
    LBP[UTR5F]->Update(log(1.0-Start[0][i]));

    LBP[UTR5F]->BestUsable(i,SwitchAny,0,&BestU);
    // Un test tordu pour casser le cou aux NaN
    if (isnan(maxi) || (BestU > maxi)) {
      maxi = BestU;
      best = -1;
    }

    // On vient de l'intergenique. 
    //
    // Deux possibilites:
    //  - on vient d'une zone entamee sur une 5' reverse (divergent genes)
    //  - on vient d'une zone entamee sur une 3' direct (flowing genes)
    //
    // On prend le meilleur des deux.

    // Sur 5' reverse
    BestU = PBest[23]-TransStartP;
    // Un test tordu pour casser le cou aux NaN
    if (isnan(maxi) || (BestU > maxi)) {
      maxi = BestU;
      best = 23;
      source = InterGen5;
      Switch = SwitchTransStart;
    }

    // Sur 3' direct
    BestU = PBest[18]-TransStartP;
    // Un test tordu pour casser le cou aux NaN
    if (isnan(maxi) || (BestU > maxi)) {
      maxi = BestU;
      best = 18;
      source = InterGen3;
      Switch = SwitchTransStart;
    }

    if (best != -1) LBP[UTR5F]->InsertNew(source,Switch,i,maxi,PrevBP[best]);

    LBP[UTR5F]->Update(log(BaseScore[8][i]));

    // ----------------------------------------------------------------
    // ---------------------- UTR 3' direct ---------------------------
    // ----------------------------------------------------------------
    maxi = NINFINITY;
    
    // On reste 3' direct
    LBP[UTR3F]->BestUsable(i,SwitchAny,0,&BestU);
    // Un test tordu pour casser le cou aux NaN
    if (isnan(maxi) || (BestU > maxi)) {
      maxi = BestU;
      best = -1;
    }
    // Ca vient d'un exon direct + STOP
    for (k = 0; k < 3; k++) {
      if ((i % 3 == k) && Stop[0][i]) {
	BestU = PBest[k+12]-StopP;
	// Un test tordu pour casser le cou aux NaN
	if (isnan(maxi) || (BestU > maxi)) {
	  maxi = BestU;
	  best = k+12;
	  Switch = SwitchStop;
	}
      }
    }
    
    if (best != -1) LBP[UTR3F]->InsertNew(best %12,Switch,i,maxi,PrevBP[best]);

    LBP[UTR3F]->Update(log(BaseScore[9][i]));

    // ----------------------------------------------------------------
    // ----------------------- UTR 5'reverse --------------------------
    // ----------------------------------------------------------------
    maxi = NINFINITY;
    
    // On reste 5' reverse
    LBP[UTR5R]->Update(log(1.0-Start[1][i]));


    LBP[UTR5R]->BestUsable(i,SwitchAny,0,&BestU);
    // Un test tordu pour casser le cou aux NaN
    if (isnan(maxi) || (BestU > maxi)) {
      maxi = BestU;
      best = -1;
    }

    // Ca vient d'un exon reverse + START
    for (k = 3; k < 6; k++) {
      if (((Data_Len-i) % 3 == k-3) && Start[1][i] != 0.0) {
	BestU = PBest[k+12]+log(Start[1][i]);
	// Un test tordu pour casser le cou aux NaN
	if (isnan(maxi) || (BestU > maxi)) {
	  maxi = BestU;
	  best = k+12;
	  Switch = SwitchStart;
	}
      }
    }

    if (best != -1) LBP[UTR5R]->InsertNew(best % 12,Switch,i,maxi,PrevBP[best]);

    LBP[UTR5R]->Update(log(BaseScore[8][i]));

    // ----------------------------------------------------------------
    // ----------------------- UTR 3'reverse --------------------------
    // ----------------------------------------------------------------
    maxi = NINFINITY;
 
    // On reste 3' reverse
    LBP[UTR3R]->BestUsable(i,SwitchAny,0,&BestU);
    // Un test tordu pour casser le cou aux NaN
    if (isnan(maxi) || (BestU > maxi)) {
      maxi = BestU;
      best = -1;
    }

    // On demarre depuis l'intergenique
    //
    // Deux possibilites:
    //  - on vient d'une zone entamee sur une 5' reverse (flowing genes)
    //  - on vient d'une zone entamee sur une 3' direct (convergent genes)
    //
    // On prend le meilleur des deux.

    // Sur 5' reverse
    BestU = PBest[24]-TransStopP;
    // Un test tordu pour casser le cou aux NaN
    if (isnan(maxi) || (BestU > maxi)) {
      maxi = BestU;
      best = 24;
      source = 12;
      Switch = SwitchTransStop;
    }

    // Sur 3' direct
    BestU = PBest[25]-TransStopP;
    // Un test tordu pour casser le cou aux NaN
    if (isnan(maxi) || (BestU > maxi)) {
      maxi = BestU;
      best = 25;
      source = 12;
      Switch = SwitchTransStop;
    }
         
    if (best != -1) LBP[UTR3R]->InsertNew(source,Switch,i,maxi,PrevBP[best]);
      
    LBP[UTR3R]->Update(log(BaseScore[10][i]));

    // ----------------------------------------------------------------
    // ---------------- Introns de phase k forward --------------------
    // ----------------------------------------------------------------
    for (k = 0; k<3; k++) {

      maxi = NINFINITY;

      // On reste intronique
      LBP[6+k]->Update(log(1.0-Acc[0][i]));
      LBP[6+k]->BestUsable(i,SwitchAny,0,&BestU);
      // Un test tordu pour casser le cou aux NaN
      if (isnan(maxi) || (BestU > maxi)) {
	maxi = BestU;
	best = -1;
      }
      // - on quitte un exon
      BestU = PBest[((i-k+3) % 3)]+log(Don[0][i]);
      // Un test tordu pour casser le cou aux NaN
      if (isnan(maxi) || (BestU > maxi)) {
	maxi = BestU;
	best = ((i-k+3) % 3);
	Switch = SwitchDon;
      }
 
      if (best != -1) LBP[6+k]->InsertNew(best,Switch,i,maxi,PrevBP[best]);
      
      LBP[6+k]->Update(log(BaseScore[6][i]));
      LBP[6+k]->Update(((ESTMatch[i] & Hit ) != 0)*EstP);

      if ((ForcedIG != NULL) && ForcedIG[i])
	LBP[6+k]->Update(IGPenalty);
    }
      
    // ----------------------------------------------------------------
    // ----------------- Introns de phase -k reverse ------------------
    // ----------------------------------------------------------------
    for (k = 0; k<3; k++) {
      
      maxi = NINFINITY;
      
      // On reste intronique
      LBP[9+k]->Update(log(1.0-Acc[1][i]));
      LBP[9+k]->BestUsable(i,SwitchAny,0,&BestU);
      // Un test tordu pour casser le cou aux NaN
      if (isnan(maxi) || (BestU > maxi)) {
	maxi = BestU;
	best = -1;
      }

      // On quitte un exon
      BestU = PBest[3+((Data_Len-i-k) % 3)]+log(Acc[1][i]);
      // Un test tordu pour casser le cou aux NaN
      if (isnan(maxi) || (BestU > maxi)) {
	maxi = BestU;
	best = 3+((Data_Len-i-k) % 3);
	Switch = SwitchAcc;
      }
      
      if (best != -1) LBP[9+k]->InsertNew(best,Switch,i,maxi,PrevBP[best]);
      
      LBP[9+k]->Update(log(BaseScore[7][i]));
      LBP[9+k]->Update(((ESTMatch[i] & Hit) != 0)*EstP);

      if ((ForcedIG != NULL) && ForcedIG[i])
	LBP[9+k]->Update(IGPenalty);
    }
  }

  LBP[ExonF1]->Update(log(ExonPrior/6.0)/2.0);
  LBP[ExonF2]->Update(log(ExonPrior/6.0)/2.0);
  LBP[ExonF3]->Update(log(ExonPrior/6.0)/2.0);
  LBP[ExonR1]->Update(log(ExonPrior/6.0)/2.0);
  LBP[ExonR2]->Update(log(ExonPrior/6.0)/2.0);
  LBP[ExonR3]->Update(log(ExonPrior/6.0)/2.0);

  LBP[IntronF1]->Update(log(IntronPrior/6.0)/2.0);
  LBP[IntronF2]->Update(log(IntronPrior/6.0)/2.0);
  LBP[IntronF3]->Update(log(IntronPrior/6.0)/2.0);
  LBP[IntronR1]->Update(log(IntronPrior/6.0)/2.0);
  LBP[IntronR2]->Update(log(IntronPrior/6.0)/2.0);
  LBP[IntronR3]->Update(log(IntronPrior/6.0)/2.0);


  // Intergenique 
  LBP[InterGen5]->Update(log(InterPrior)/2.0); 
  LBP[InterGen3]->Update(log(InterPrior)/2.0); 

  // UTR 5' et 3'
  LBP[UTR5F]->Update(log(InterPrior)/2.0);
  LBP[UTR3F]->Update(log(InterPrior)/2.0);  
  LBP[UTR5R]->Update(log(InterPrior)/2.0);
  LBP[UTR3R]->Update(log(InterPrior)/2.0);

  for (i = 0; i < 18; i++) {
    PrevBP[i] = LBP[i]->BestUsable(Data_Len+1,SwitchAny,0,&BestU);    
    LBP[i]->InsertNew(i,0xFF,Data_Len+1,BestU,PrevBP[i]);
  }
  
  // Backtrace in the graph
  LBP[0]->BestUsable(Data_Len+1,SwitchAny,0,&maxi);
  j = 0;

  for (i = 1; i < 18 ; i++) {
    LBP[i]->BestUsable(Data_Len+1,SwitchAny,0,&BestU);
    // Un test tordu pour casser le cou aux NaN
    if (isnan(maxi) || (BestU > maxi)) {
      maxi = BestU;
      j = i;
    }
  }

  LBP[j]->BackTrace(Choice);
  
    for  (i = 0;  i < 18;  i ++) LBP[i]->Zap();

  if (!PorteOuverte && Data_Len > 6000) 
    exit(2);
  
  fprintf(stderr,"Optimal path length = %#f\n",maxi);

  if (graph) PlotPredictions(Data_Len,Choice,Stop,Start,Acc,Don);

#include "Output.h"

  // free used memory
  if (graph) ClosePNG();

  delete TheSeq;
  
  delete [] Stop[0];
  delete [] Stop[1];

  delete [] ESTMatch;
  
  delete [] Start[0];
  delete [] Start[1];

  delete [] Don[0];
  delete [] Don[1];
  
  delete [] Acc[0];
  delete [] Acc[1];

   delete [] Choice;

  for  (i = 0;  i < 11;  i ++) delete [] BaseScore[i];
  }

  // free remaining used memory
  if (IMMatrix[5] != IMMatrix[4]) delete  IMMatrix[5];

  for  (i = 0;  i < 5;  i ++)
    delete  IMMatrix[i];
    
  return  0;
}
