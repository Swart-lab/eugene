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
// Fri Mar 29 22:17:29 2002: kludge: les penalites de start sont ignorees des UTR
//                             si l'on a est dans un intron UTR d'apres les hits EST

// Mon May 13 09:52:35 2002: PAYFORIGNORE flag: indique si on paye
//                              quand on ne prend pas un aiguillage.
// TODO:
// supprimer Choice
// remettre Frameshifts
// alignement cds cDNA et proteines sur les splice sites

// supprimer l'approximation longueur single (on peut rater un bon
// epissage si un meilleur START l'occulte mais qu'il est trop pres
// pour etre utilisable). Il faudrait faire 3 pistes single + 3 pistes
// non single.

#define  VERSION "1.2b (280802)"
#define  VERSION_PAR "21_08_02"
#define PAYTOIGNORE 1

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <vector>
#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cctype>
#include <climits>
#include <cfloat>
#include <ctime>
#include <cassert>
#include <cerrno>
#include <unistd.h>

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
#include "structure.h"
#include "yacc.tab.h"

//#include "dmalloc.h"
//#include "dmallocc.cc"
//#define STAND 1

#define DFT_MATRIX          "default.mat"
#define DFT_OUTPUT          "./"

#ifndef FILENAME_MAX
#define FILENAME_MAX        1024
#endif



// ------------------ Globals ---------------------
// Les etats
  
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

// Les Hits EST

const unsigned char HitForward    = 0x1;
const unsigned char MLeftForward  = 0x2;
const unsigned char GapForward    = 0x4;
const unsigned char MRightForward = 0x8; 

// shift to go from Hits to ...
const  unsigned int HitToMLeft  = 1;
const  unsigned int HitToGap    = 2;
const  unsigned int HitToMRight = 3;

const unsigned char HitReverse    = 0x10;
const unsigned char MLeftReverse  = 0x20;
const unsigned char GapReverse    = 0x40;
const unsigned char MRightReverse = 0x80; 

const unsigned char Hit           = HitForward    | HitReverse;
const unsigned char MLeft         = MLeftForward  | MLeftReverse;
const unsigned char MForward      = MLeftForward  | MRightForward;
const unsigned char MReverse      = MLeftReverse  | MRightReverse;
const unsigned char Gap           = GapForward    | GapReverse;
const unsigned char MRight        = MRightForward | MRightReverse;
const unsigned char Margin        = MRight        | MLeft;

const unsigned char NotAHit       = Margin | Gap;

#define Inconsistent(x) (((x) & Hit) && ((x) & Gap))

double Score [9],NScore[9];
int Data_Len;
int normopt,blastopt,estopt,estanal,ncopt,raflopt,userinfo;
int window, offset,graph,resx,resy;
int    gfrom,gfromSave,gto,gtoSave,golap,glen;

extern char     *optarg;   
extern int       optind;

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

// RAFL: Riken Arabidopsis Full Length cDNA
typedef struct RAFLgene{
  int deb;
  int fin;
  signed char sens;
  char ID[FILENAME_MAX+1];
} RAFLgene;

inline bool Before(const RAFLgene &A,const RAFLgene &B) {
  return (A.deb < B.deb);
}

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


//--------------------------------------------------
int HitsCompareLex(const void *A, const void *B)
{
  Hits **UA,**UB;
  
  UA = (Hits **) A;
  UB = (Hits **) B;

  if ((*UA)->Start < (*UB)->Start) return -1;
  if ((*UA)->Start > (*UB)->Start) return 1;
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

//--------------------------------------------------
// Convertit les phases 0-6 en 1 2 3 -1 -2 -3 0
//--------------------------------------------------
int PhaseAdapt(char p)
{
  if (p >= 12) return 0;
  else if (p < 3) return (1+p);
  else if (p < 6) return (2-p);
  else if (p < 9) return (p-2);
  else return (5-p);
  
}

//--------------------------------------------------
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

  case 17:
    printf("IG ");
    return;

  default:
    fprintf(stderr,"ERROR: unexpected Choice value\n");
    exit(1);
  }
  return;
}

//--------------------------------------------------
// Convertit les phases 1 2 3 -1 -2 -3 0 en 0-6
//--------------------------------------------------
char ph06(char p)
{
  if (p == 0) return 6;
  else if (p > 0) return (p-1);
  else return 2-p;   
}
//--------------------------------------------------
// Dump les signaux au format util user
//--------------------------------------------------
void DumpSignals(int Len, REAL** Donor, REAL **Acceptor,REAL** Stop, REAL **Start,FILE* flot)
{
  int i;
  for (i=0; i< Len; i++) {
    //    if (Stop[0][i]) fprintf(flot,"stop f %d %a\n",i+1,Stop[0][i]);
    //    if (Stop[1][i]) fprintf(flot,"stop r %d %a\n",i+1,Stop[1][i]);

    if (Start[0][i]) fprintf(flot,"start f %d %a\n",i+1,Start[0][i]);
    if (Start[1][i]) fprintf(flot,"start r %d %a\n",i+1,Start[1][i]);

    if (Donor[0][i]) fprintf(flot,"donor f %d %a\n",i+1,Donor[0][i]);
    if (Donor[1][i]) fprintf(flot,"donor r %d %a\n",i+1,Donor[1][i]);

    if (Acceptor[0][i]) fprintf(flot,"acceptor f %d %a\n",i+1,Acceptor[0][i]);
    if (Acceptor[1][i]) fprintf(flot,"acceptor r %d %a\n",i+1,Acceptor[1][i]);
  }
}
//--------------------------------------------------
Hits **ESTAnalyzer(FILE *ESTFile,   unsigned char    *ESTMatch,
		 int EstM, REAL** Don, REAL **Acc, int *NumEST)
{
  int i,j,k;
  int Rejected = 0;
  Hits *ThisEST = NULL,*AllEST = NULL;
  Block *ThisBlock = NULL;
  
  fprintf(stderr,"Reading cDNA hits...");
  fflush(stderr);

  AllEST = AllEST->ReadFromFile(ESTFile,NumEST);
  
  fprintf(stderr,"%d sequences read.\n",*NumEST);
  fflush(stderr);

  // on trie les hits sur le nombre de gaps et la
  // longueur. L'idee est d'eliminer les epissages partiels et de
  // favoriser la longueur de toute facon.
  
  Hits **HitTable = new Hits *[*NumEST+1];
  
  for (i = 0, ThisEST = AllEST; i < *NumEST; i++, ThisEST = ThisEST->Next)
    HitTable[i] = ThisEST;

  // pour memorise le premier (a liberer)
  HitTable[*NumEST] = AllEST;
  
  qsort((void *)HitTable,*NumEST,sizeof(void *),HitsCompare);
  
  for (int index = 0; index < *NumEST; index++) {
    
    int Inc;
    int ExonInc,TheStrand;
    REAL WorstSpliceF, WorstSpliceR;
    REAL DonF,AccF,DonR,AccR;
    
    ThisEST = HitTable[index];

    // Le veritable  brin est a priori indetermine
    TheStrand = HitForward | HitReverse;
    Inc = 0;
    ExonInc = 0;
    WorstSpliceF = WorstSpliceR = 1.0;
    
    
    // Look for each match in the current Hit
    ThisBlock = ThisEST->Match;

    // First Step: tries to determine strand
    while (ThisBlock) {
      // si ona un gap ?
      if ((ThisBlock->Prev != NULL) && 
	  abs(ThisBlock->Prev->LEnd - ThisBlock->LStart) <= 6) {
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

	//	printf("Extreme splices: %f %f - %f %f\n",DonF,AccF,DonR,AccR);
	WorstSpliceF = Min(WorstSpliceF,DonF);
	WorstSpliceF = Min(WorstSpliceF,AccF);
	WorstSpliceR = Min(WorstSpliceR,DonR);
	WorstSpliceR = Min(WorstSpliceR,AccR);
      }
      
      ThisBlock = ThisBlock->Next;
      }
    
    // Tous les blocs ont ete traites
    if (WorstSpliceF == 0.0) TheStrand &= (~HitForward);
    if (WorstSpliceR == 0.0) TheStrand &= (~HitReverse);

    //    printf("Strand %d\n",TheStrand);

    // next iteration on the same EST
    ThisBlock = ThisEST->Match;
    while (TheStrand && ThisBlock) {
      // Check for consistency with already read Hits
      // The Inc flag will keep the Inconsistency status
      // 1 - inconsistent with a previous EST
      // 2 - no splice site on the borders of a gap
      // 3 - an exon contains STRONG donor on both strands

      for (i = ThisBlock->Start+EstM; !Inc && i <= ThisBlock->End-EstM; i++) {

	if (((ESTMatch[i] & Hit) && !(ESTMatch[i] & TheStrand)) |
	    (Inconsistent(ESTMatch[i] | TheStrand))) {
	  fprintf(stderr,"   [%s]: inconsistent hit [%d-%d]\n",
		  ThisEST->Name,ThisBlock->Start+1,ThisBlock->End+1);
	  Inc = 1;
	}
      }
      
      // si on a un gap
      if ((ThisBlock->Prev != NULL) && abs(ThisBlock->Prev->LEnd - ThisBlock->LStart) <= 6) {
        
        for (i = ThisBlock->Prev->End+1+EstM; !Inc && i < ThisBlock->Start-EstM; i++) 
	  if (((ESTMatch[i] & Gap) && !(ESTMatch[i] & (TheStrand << HitToGap))) |
	    (Inconsistent(ESTMatch[i] | (TheStrand << HitToGap)))) {
            fprintf(stderr,"   [%s]: inconsistent gap [%d-%d]\n",
                    ThisEST->Name,ThisBlock->Prev->End+2,ThisBlock->Start);
            Inc = 1;
          }
      }

      DonF = 0.0;
      DonR = 0.0;
      
      // calcul des sites d'epissage internes a l'exon
      for (i = ThisBlock->Start+EstM+1; !Inc && i <= ThisBlock->End-EstM-1; i++) {
	DonF = Max(DonF,Don[0][i]);
	DonR = Max(DonR,Don[1][i]);
	//	Don[0][i] *= 0.99;
	//	Don[1][i] *= 0.99;
      }
      //      printf("Internal splices: %f - %f\n",DonF,DonR);
      if (DonF > 0.95) ExonInc |= 1;
      if (DonR > 0.95) ExonInc |= 2;
      if (ExonInc == 3 && !Inc && !ThisEST->NGaps) {
	fprintf(stderr,"   [%s]: Gapless EST with strong donor [%d-%d]\n",
		    ThisEST->Name,ThisBlock->Start+1,ThisBlock->End+1);
	Inc = 3;
      }
      
      ThisBlock = ThisBlock->Next;
    }
    
    if (!TheStrand) {
      fprintf(stderr, "   [%s]: no matching splice site\n",ThisEST->Name);
      Inc =2;
    }
    
    // Si une incoherence est detectee, on va jeter la sequence
    
    if (Inc) {
      Rejected++;
      ThisEST->Rejected = TRUE;
      
      ThisBlock = ThisEST->Match;
      while (ThisBlock) {
	for (i = ThisBlock->Start; i<= ThisBlock->End; i++) {
	  if (TheStrand & HitForward) PlotBarI(i, 4,0.9,1,7);
	  if (TheStrand & HitReverse) PlotBarI(i,-4,0.9,1,7);
	}
	
	if (graph && (ThisBlock->Prev != NULL) && 
	    abs(ThisBlock->Prev->LEnd-ThisBlock->LStart) <= 6) {
	  if (TheStrand & HitForward) 
	    PlotLine(ThisBlock->Prev->End,ThisBlock->Start,4,4,0.9,0.9,7);
	  if (TheStrand & HitReverse) 
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
	
	for (i = LBoundary; i <= RBoundary; i++) {
	  if (ESTMatch[i] & Hit)
	    ESTMatch[i] |= (ESTMatch[i] & TheStrand);
	  else
	    ESTMatch[i] |= TheStrand;
	}

	if (ThisBlock->Prev == NULL)
	  for (i = ThisBlock->Start; i<LBoundary; i++)
	    ESTMatch[i] |= TheStrand << HitToMLeft;

	if (ThisBlock->Next == NULL)
	  for (i = RBoundary+1; i <= ThisBlock->End; i++)
	    ESTMatch[i] |= TheStrand << HitToMRight;
	
	if ((ThisBlock->Prev != NULL) && 
	    abs(ThisBlock->Prev->LEnd-ThisBlock->LStart) <= 6) {
	  
	  for (i = Max(0,ThisBlock->Prev->End+1-EstM); 
	       i < Min(Data_Len, ThisBlock->Prev->End+1+EstM); i++) 
	    if (ESTMatch[i] & MLeft)
	      ESTMatch[i] |= (ESTMatch[i] & (TheStrand << HitToMLeft));
	    else
	      ESTMatch[i] |= TheStrand << HitToMLeft;
	  
	  for (i = ThisBlock->Prev->End+1; i < ThisBlock->Start-1; i++)
	    if (ESTMatch[i] & Gap)
	      ESTMatch[i] |= (ESTMatch[i] & (TheStrand << HitToGap));
	    else
	    ESTMatch[i] |= TheStrand << HitToGap;
	  
	  for (i =  Max(0,ThisBlock->Start-1-EstM);
	       i <  Min(Data_Len, ThisBlock->Start-1+EstM); i++)
	    if (ESTMatch[i] & MRight)
	      ESTMatch[i] |= (ESTMatch[i] & (TheStrand << HitToMRight));
	    else
	      ESTMatch[i] |= TheStrand << HitToMRight;

	  
	  if (graph) {
	    if (TheStrand & HitForward) 
	      PlotLine(ThisBlock->Prev->End,ThisBlock->Start,4,4,0.6,0.6,2);
	    if (TheStrand & HitReverse) 
	      PlotLine(ThisBlock->Prev->End,ThisBlock->Start,-4,-4,0.6,0.6,2);
	  }
	  
	}
	ThisBlock = ThisBlock->Next;
      }
    }
  }
  
  if (graph)
    for (i = 0; i < Data_Len; i++) {
      if (ESTMatch[i] & HitForward) 
	PlotBarI(i, 4,0.6,1,2);
      if (ESTMatch[i] & HitReverse) 
	PlotBarI(i, -4,0.6,1,2);
    }

  if (Rejected)
    fprintf(stderr,"A total of %d/%d sequences rejected\n",Rejected,*NumEST);
  return HitTable;
}

// -------------------------------------------------------------------------
// Verif coherence EST: calcul le nombre de nuc. coherents et
// incoherents avec les match est
// debut/fin/etat: debut et fin de la seq. dont l'etat est etat
// Match: resume des hits EST
// cons/incons: retour des valeurs
// -------------------------------------------------------------------------
void CheckConsistency(int debut, int fin, int etat, 
		      unsigned char *Match, int * cons, int* incons)
{
  int i, con = 0, inc = 0;

  // les valeurs qui sont coherentes avec chaque etat
  const unsigned char Consistent[18] = {
    HitForward|MForward,    HitForward|MForward,    HitForward|MForward,
    HitReverse|MReverse,    HitReverse|MReverse,    HitReverse|MReverse,
    GapForward|MForward,    GapForward|MForward,    GapForward|MForward,
    GapReverse|MReverse,    GapReverse|MReverse,    GapReverse|MReverse,
    0,
    HitForward|MForward|GapForward,HitForward|MForward|GapForward,
    HitReverse|MReverse|GapReverse,    HitReverse|MReverse|GapReverse,
    0
  };
  
 const unsigned char MaskConsistent[18] = {
    Hit|Margin,    Hit|Margin,    Hit|Margin,    
    Hit|Margin,    Hit|Margin,    Hit|Margin,    
    Gap|Margin,    Gap|Margin,    Gap|Margin,
    Gap|Margin,    Gap|Margin,    Gap|Margin,
    0,
    Margin|Hit|Gap,    Margin|Hit|Gap,
    Margin|Hit|Gap,    Margin|Hit|Gap,
    0
  };

 if (debut == -1) debut = 0;

 for (i = debut; i <fin; i++) {
   // y a t'il de l'info
   if (Match[i]) {
     // y a t'il une info incoherente avec l'etat
     if (Match[i] & ~MaskConsistent[etat]) 
       inc++;
     else if (Match[i] & Consistent[etat]) 
       con++;
     else 
       inc++;
   }
 }

 *cons = con;
 *incons = inc;
}
// -------------------------------------------------------------------------
// Tdebut/fin = debut/fin de transcript
// debut/fin = debut/fin de traduit
// -------------------------------------------------------------------------
void ESTSupport(char * Choice, int Tdebut, int Tfin, int debut,int fin,  Hits **HitTable, int Size)
{
  static int EstIndex;
  int supported = 0;
  int CDSsupported = 0;
  unsigned char *Sup;
  Block *ThisBlock;
  int ConsistentEST,i;
  int from = 0, to = 0, ESTEnd = 0;
  
  if (Choice == NULL) { 
    EstIndex = 0;
    return;
  }

  Sup = new unsigned char[Tfin-Tdebut+1]; 

  for (i=0; i <= Tfin-Tdebut; i++)
    Sup[i]=0;
  
  // si la fin du codant n'a pas ete rencontree
  if (fin == -1) fin = Tfin;
  if ((debut == -1) || (debut > Tfin)) debut = Tfin+1;

  while (EstIndex < Size) {
    // le dernier transcrit exploitable est passe
    if (HitTable[EstIndex]->Start > Tfin) break;

    ConsistentEST = 1;
    ThisBlock = HitTable[EstIndex]->Match;

     while (ThisBlock && ConsistentEST) {
       // si on a un gap
       if (ThisBlock->Prev != NULL) {
	 // intersection
	 from = Max(Tdebut,ThisBlock->Prev->End+1);
	 to = Min(Tfin,ThisBlock->Start-1);
	 
	 for (i = from; i <=to; i++)
	   if ((Choice[i+1] < IntronF1) || Choice[i+1] == InterGen5 || Choice[i+1] == InterGen3)
	     ConsistentEST = 0;
       }

       from = Max(Tdebut,ThisBlock->Start);
       to = Min(Tfin,ThisBlock->End);
       ESTEnd = ThisBlock->End;


       for (i = from; i <= to; i++)
	 if (((Choice[i+1] > ExonR3) && (Choice[i+1] <= InterGen5)) || Choice[i+1] == InterGen3)
	   ConsistentEST = 0;

       ThisBlock = ThisBlock->Next;
     }
      

     printf("cDNA  %-12s %7d %7d     %4d     %2d introns    ",HitTable[EstIndex]->Name,
	    HitTable[EstIndex]->Start+1,ESTEnd+1,
	    HitTable[EstIndex]->Length,HitTable[EstIndex]->NGaps);

     if (HitTable[EstIndex]->Rejected) printf("Filtered ");
     else if (!ConsistentEST) printf("Inconsistent");
     else {
       if ((HitTable[EstIndex]->Start) <= Tdebut && (ESTEnd >= Tfin))
	 printf("Full Transcript Support");
       else if ((HitTable[EstIndex]->Start) <= debut && (ESTEnd >= fin))
	 printf("Full Coding Support");
       else printf("Support");
       

       ThisBlock = HitTable[EstIndex]->Match;
       while (ThisBlock) {
	 if (ThisBlock->Prev != NULL) {
	   // intersection
	   from = Max(Tdebut,ThisBlock->Prev->End+1);
	   to = Min(Tfin,ThisBlock->Start-1);

	   for (i= from; i <= to; i++) 
	     if (!Sup[i-Tdebut]) {
	       Sup[i-Tdebut] = TRUE;
	       supported++;
	       if ((i >=debut) && (i <=fin)) CDSsupported++;
	     }
	 }

	 from = Max(Tdebut,ThisBlock->Start);
	 to = Min(Tfin,ThisBlock->End);

	 for (i= from; i <= to; i++) 
	   if (!Sup[i-Tdebut]) {
	     Sup[i-Tdebut] = TRUE;
	     supported++;
	     if ((i >=debut) && (i <=fin)) CDSsupported++;
	   }
	 ThisBlock = ThisBlock->Next;
       }
     }
     printf("\n");
     EstIndex++;
  }
  if (fin >= debut)
  printf("      CDS          %7d %7d    %5d     supported on %d bases\n",
	 debut+1,fin+1,fin-debut+1,CDSsupported);
  printf("      Gene         %7d %7d    %5d     supported on %d bases\n",
	 Tdebut+1,Tfin+1,Tfin-Tdebut+1,supported);
  delete Sup;
  return;
}

// -------------------------------------------------------------------------
//            MAIN
// -------------------------------------------------------------------------
int main  (int argc, char * argv [])
{
  int              i, j, k, carg, errflag;
  FILE             *fp;
  REAL             BaseScore[13],*Start[2] ;
  REAL             *Don[2],*Acc[2];
  REAL             *Stop[2];
  unsigned char    *ESTMatch = NULL;

  REAL             *ProtMatch = NULL;
  REAL             *ProtMatchLevel = NULL;
  int              *ProtMatchPhase = NULL;

  BString_Array    *IMMatrix[7];
  char             clef[20];
  DNASeq           *TheSeq;
  char             printopt;
  char             outputname[FILENAME_MAX+1];
  char             parname[FILENAME_MAX+1];
  char             blastArg[FILENAME_MAX+1];
  char             matname[FILENAME_MAX+1], tempname[FILENAME_MAX+1], fstname[FILENAME_MAX+1];
  char             grname[FILENAME_MAX+1];
  char             grnameArg[FILENAME_MAX+1]; // argument -g"...."
  char             *grnameFile = NULL;
  int              EstM;
  char             versionPAR[FILENAME_MAX+1];
  double           FsP,StartP,StartB,StopP,TransStartP,TransStopP;
  double           AccP[2],AccB[2],DonP[2],DonB[2],BlastS[8],EstP;
  double           TransitionWeight[5];
  
  // Nombre de paramètre dans le fichier .par
  const int NBPARAM  = 39;
  const int MAX_LINE = 300;
  char      Line [MAX_LINE];

  unsigned char    *ForcedIG  = NULL;
  Hits             **HitTable = NULL;
  int              NumEST;
  
  char             *EugDir;

  std::vector <RAFLgene> RAFL;
  int              RAFLpos = 0;  //position par rapport a un gene RAFL
  int              RAFLindex = 0;   // index du RIKEN en cours
  const REAL RAFLPenalty = NINFINITY; 
  
  // prior on the initial state, selon Sato et al 1999 / Terryn et
  // al. 1999
  REAL ExonPrior = 0.33,IntronPrior = 0.17;
  REAL InterPrior = 0.4,FivePrimePrior = 0.03,ThreePrimePrior = 0.07;

  const REAL DontCrossStop = NINFINITY;
  const REAL IGPenalty = -1.0; 

  // Les longueurs. 

  int MinFivePrime, MinThreePrime;
  int MinEx, MinIn, MinSg,MinFlow,MinConv,MinDiv;
  int MinLength[18];

  // process args 
  // default values

  glen = -1;
  golap = -1;
  gfrom = -1;
  gfromSave = -1;
  gto = -1;
  gtoSave = -1;
  resx = 900;                 // x res for graph. output
  resy = 400;                 // y res for graph. output
  graph = FALSE;              // don't produce a graphical output
  estopt = FALSE;             // don't try to read a EST file
  estanal = FALSE;            // don't try to analyze EST support
  userinfo = FALSE;           // shall we read a user info file
  raflopt = FALSE;            // don't try to read a est.rafl file
  blastopt = FALSE;           // don't try to read a blast file
  ncopt = FALSE;              // don't try to read a non coding input
  normopt = 1;                // normalize across frames
  window  = 97;               // window length
  printopt = 'l';             // short print format 
  offset = 0;                 // no offset

  (void) strcpy(matname, DFT_MATRIX); // default matrix file
  (void) strcpy(outputname, DFT_OUTPUT); // default output
  *fstname = '\000';                  // no default input
  parname[0] = '0';                   // no -P
  errflag = 0;

  while ((carg = getopt(argc, argv, "UREdrshm:w:f:n:o:p:x:y:c:u:v:g::b::l:P:O:")) != -1) {
    
    switch (carg) {
      
    case 'n':            /* -n normalize across frames      */
      if (! GetIArg(optarg, &normopt, normopt))
	errflag++;
      break;
      
    case 'h':           /* help                             */
      errflag++;
      break;

    case 's':           /* Single gene mode. Start/End on IG only   */
      ExonPrior = 0.0;
      IntronPrior = 0.0;
      //      InterPrior = 1.0;
      FivePrimePrior = 0.0;
      ThreePrimePrior = 0.0;
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
      gfromSave = gfrom;
      break;
      
    case 'v':            /* -v To      */
      if (! GetIArg(optarg, &gto, gto))
	errflag++;
      gtoSave = gto;
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
      if (optarg) {
	if (strspn(".fasta", optarg)    != 6 
	    && strspn(".fsa", optarg)   != 4
	    && strspn(".tfa", optarg)   != 4
	    && strspn(".txt", optarg)   != 4
	    && optarg[0] != '-') {
	  if(strcmp(optarg, "NO"))      // Argument suivant le -g
	    strcpy(grnameArg,optarg);   // != NO pris en compte
	  else 
	    grnameArg[0] = 0;           // == NO non pris en compte
	  graph = TRUE;
	}
	else errflag++;
      }
      else
	grnameArg[0] = 0;
      break;
      
    case 'p':           /* print opt: short/long/detailed   */
      printopt = optarg[0];
      if ((printopt == 'h') &&  graph == FALSE) {
	graph = TRUE; // HTML output means graphical output 
	grnameArg[0] = 0;
      }
      if ((printopt != 's') && (printopt != 'l') && (printopt != 'g') && 
	  (printopt != 'd') && (printopt != 'h'))
	errflag++;
      break;
      
    case 'm':           /* -m matrix                        */
      (void) strcpy(matname, optarg);
      break;

    case 'O':           /* -O output                        */
      (void) strcpy(outputname, optarg);
      break;

    case 'P':           /* -P .par                          */
      (void) strcpy(parname, optarg);
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
      if (optarg) {
	(void) strcpy(blastArg, optarg);
	blastopt = TRUE;
	if(strlen(blastArg) > 8)
	  errflag++;
      }
      else {
	blastArg[0] = '1';
	blastArg[1] = '2';
	blastArg[2] = '3';
      }
      break;
      
    case 'E':
      estanal = TRUE;
      break;
      
    case 'U':
      userinfo = TRUE;
      break;

    case 'd':           /* -d  use cDNA blastn results      */
      estopt = TRUE;
      break; 

    case 'R':           /* -R use RAFL-like EST*/
      raflopt = TRUE;
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
    fprintf(stderr, "Usage: EuGene [-h] [-m matrix] [-P .par] [-n 0|1|2] [-s] [-p h|g|s|l|d]\n"
	    "              [-w window] [-b {levels}] [-d] [-R] [-E] [-U] [-o offset]\n"
	    "              [-g {graphArg}] [-u start] [-v end] [-l len] [-c olap]\n"
	    "              [-x xres] [-y yres] FASTA files\n");
    exit(1);
  }
  
  fprintf(stderr,"EuGene rel. %s\n",VERSION);
  fprintf(stderr,"Loading parameters file...");
  fflush(stderr);

  // Lecture des parametres (.par)
  EugDir = getenv("EUGENEDIR");
  
  if(parname[0] == '0') {
    strcpy(tempname,argv[0]);
    strcat(tempname,".par");
    fp = FileOpen(EugDir,BaseName(tempname),"r");
  }
  else {
    fp = FileOpen(EugDir,parname,"r");
    strcpy(tempname,parname);
  }
  
  i = 0;
  fgets (Line, MAX_LINE, fp);
  fgets (Line, MAX_LINE, fp);
  fgets (Line, MAX_LINE, fp);
  if (fgets (Line, MAX_LINE, fp) != NULL)
  if (sscanf(Line, "%s", clef))
    if (fgets (Line, MAX_LINE, fp) != NULL)
    if (sscanf(Line, "%s", versionPAR))
      if (fgets (Line, MAX_LINE, fp) != NULL)
      if (sscanf(Line, "%lf", &FsP))
        if (fgets (Line, MAX_LINE, fp) != NULL)
        if (sscanf(Line, "%lf %lf", &StartP, &StartB) == 2)
          if (fgets (Line, MAX_LINE, fp) != NULL)
          if (sscanf(Line, "%lf", &StopP))
            if (fgets (Line, MAX_LINE, fp) != NULL)
            if (sscanf(Line, "%lf", &TransStartP))
              if (fgets (Line, MAX_LINE, fp) != NULL)
              if (sscanf(Line, "%lf", &TransStopP))
       	        if (fgets (Line, MAX_LINE, fp) != NULL)
                if (sscanf(Line, "%lf %lf", &AccP[0], &AccB[0]) == 2)
                  if (fgets (Line, MAX_LINE, fp) != NULL)
     	          if (sscanf(Line, "%lf %lf", &DonP[0], &DonB[0]) == 2)
                    if (fgets (Line, MAX_LINE, fp) != NULL)
    	            if (sscanf(Line, "%lf %lf", &AccP[1], &AccB[1]) == 2)
                      if (fgets (Line, MAX_LINE, fp) != NULL)
    	              if (sscanf(Line, "%lf %lf", &DonP[1], &DonB[1]) == 2)
			if (fgets (Line, MAX_LINE, fp) != NULL)
			if (sscanf(Line, "%lf", &EstP))
			  if (fgets (Line, MAX_LINE, fp) != NULL)
			  if (sscanf(Line, "%d", &EstM))
			    if (fgets (Line, MAX_LINE, fp) != NULL)
			    if (fgets (Line, MAX_LINE, fp) != NULL)
			    if (sscanf(Line, "%lf %lf %lf %lf %lf %lf %lf %lf",
				       &BlastS[0],&BlastS[1],&BlastS[2],&BlastS[3],
				       &BlastS[4],&BlastS[5],&BlastS[6],&BlastS[7]) == 8)
			      if (fgets (Line, MAX_LINE, fp) != NULL)
			      if (fgets (Line, MAX_LINE, fp) != NULL)
			      if (sscanf(Line, "%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d",
					 &MinEx, &MinIn, &MinSg, &MinFlow, &MinConv,
					 &MinDiv, &MinFivePrime, &MinThreePrime) == 8)
				if (fgets (Line, MAX_LINE, fp) != NULL)
				if (fgets (Line, MAX_LINE, fp) != NULL)
				if (sscanf(Line, "%lf %lf %lf %lf %lf",
					   &TransitionWeight[0], &TransitionWeight[1],
					   &TransitionWeight[2], &TransitionWeight[3],
					   &TransitionWeight[4]) == 5) {
				  i = NBPARAM;
				  fprintf(stderr,"done\n");
				  if(strcmp(versionPAR,VERSION_PAR)) {
				    fprintf(stderr,"Incorrect parameter file version : %s\n",
					    versionPAR);
				    fprintf(stderr,"Version %s required\n", VERSION_PAR);
				    exit(2);
				  }
				  else
				    fprintf(stderr,"Parameters file %s\n", versionPAR);
				}
  if (i != NBPARAM)
    {
      if (EugDir != NULL)
	fprintf(stderr, "Incorrect parameter file %s/%s\n",EugDir,BaseName(tempname));
      else
	fprintf(stderr, "Incorrect parameter file %s\n",BaseName(tempname));
      exit(2);
    }
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

  // ---------------------------------------------------------------------------  
  // Lecture des modeles de Markov
  // ---------------------------------------------------------------------------  

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

  // On essaie ensuite de lire un 6eme modele. Si cela echoue,
  // le modele interG est utilise pour les UTR
  IMMatrix[6] = new BString_Array(MODEL_LEN, ALPHABET_SIZE);
  if (IMMatrix[6]->Read(fp)) {
    fprintf(stderr,"- No UTR model found, using introns model. ");
    delete IMMatrix[6];
    IMMatrix[6] = IMMatrix[3];
    IMMatrix[5] = IMMatrix[3];
  } else {
    fprintf(stderr,"6 ");
    IMMatrix[5] = new BString_Array(MODEL_LEN, ALPHABET_SIZE);
    if (IMMatrix[5]->Read(fp)) {
      fprintf(stderr,"- No second UTR model found, using intron model. ");
      delete IMMatrix[5];
      IMMatrix[5] = IMMatrix[3];
    } else fprintf(stderr,"7 ");
  }
  
  fprintf(stderr,"done\n");
  fclose(fp);

  // ---------------------------------------------------------------------------  
  // Lecture de la sequence
  // ---------------------------------------------------------------------------  
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
  
  // ---------------------------------------------------------------------------  
  // Allocation match prot/est
  // ---------------------------------------------------------------------------  
  
  if (estopt) {
    ESTMatch = new unsigned char[Data_Len+1];
    for (i = 0; i<= Data_Len; i++) 
      ESTMatch[i] = 0;
  }
  
  if (blastopt) {
    ProtMatch = new REAL[Data_Len+1];
    ProtMatchLevel = new REAL[Data_Len+1];
    ProtMatchPhase = new int[Data_Len+1];
    for (i = 0; i<= Data_Len; i++) {
      ProtMatch[i] = 0.0;
      ProtMatchLevel[i] = 0.0;
      ProtMatchPhase[i]=0;
    }
  }
  // ---------------------------------------------------------------------------
  //Acquisition des signaux
  // ---------------------------------------------------------------------------
  Stop[0] = new REAL[Data_Len+1];
  Stop[1] = new REAL[Data_Len+1];
  
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

  Check_Start(TheSeq,Start);
  
  fprintf(stderr, "Reading splice site files...");  
  fflush(stderr);
  
  Read_Splice(fstname,1,Data_Len,Acc[0],Don[0],AccP,AccB,DonP,DonB);
  fprintf(stderr,"forward,");
  fflush(stderr);
  
  Read_Splice(fstname,-1,Data_Len, Acc[1], Don[1],AccP,AccB,DonP,DonB);
  fprintf(stderr," reverse done\n");

   Check_Splices(TheSeq,Acc,Don);

  // ---------------------------------------------------------------------------
  // Lecture donnees user
  // ---------------------------------------------------------------------------
  if (userinfo) {
    for (i = 0; i < 13; i++) Weights[i] = BaseScore+i;
    Weights[13] = Start[0];    Weights[14] = Start[1];
    Weights[15] = Stop[0];    Weights[16] = Stop[1];
    Weights[17] = Acc[0];    Weights[18] = Acc[1];
    Weights[19] = Don[0];    Weights[20] = Don[1];
    fprintf(stderr,"Loading user data...");
    
    strcpy(tempname,fstname);
    strcat(tempname,".user");
    errflag = Utilisateur(tempname);   //prise en compte de donnees utilisateur
    if (errflag) {
      userinfo = FALSE;
      fprintf(stderr,"none found\n");
    }
    else fprintf(stderr,"done\n");
  }
  
  if (userinfo) {
    UserInfoList = SignalUser;
    //    WriteUtils(UserInfoList,stdout);
    
    for (i = 0; i <= Data_Len; i++) {
      //application des regles utilisateur pour i
      Util(i,UserInfoList);
      for (k = 13; k<21; k++) Weights[k]++;
    }
  }
  // ---------------------------------------------------------------------------
  // Preparation sortie graphique + Scores
  // ---------------------------------------------------------------------------
  if (graph) {
    // Récupération du nom du fichier d'origine
    grnameFile = BaseName(fstname);

    // Construction du nom de sortie (*.png)
    strcpy(grname,outputname);
    strcat(grname,grnameFile);
    *rindex(grname,'.') = 0;   // on enleve l'extension (.fasta typ.)
    if(grnameArg[0] != 0) {    // -p h sans -g ou -g NO pas d'argument
      strcat(grname,".");
      strcat(grname,grnameArg);
    }
    
    if ((gfrom <= 0)|| (gfrom >= Data_Len)) gfrom = 1;
    if ((gto <= 0)  || (gto <= gfrom) || (gto > Data_Len)) gto = Data_Len;
    
    if (gfromSave != -1) {
      sprintf(grname+strlen(grname), ".%d", gfrom);
      if (gtoSave == -1) sprintf(grname+strlen(grname), "-%d", Data_Len);
    }
    if (gtoSave != -1) {
      if (gfromSave == -1) sprintf(grname+strlen(grname), ".%d", 1);
      sprintf(grname+strlen(grname), "-%d", gto);
    }
    
    gfrom--;
    gto--;
    
    if (glen < 0) glen = ((gto-gfrom+1 <= 6000) ? gto-gfrom+1 : 6000);
    
    InitPNG(resx,resy,offset,gfrom,gto,golap,glen,grname);
    PlotScore(TheSeq,IMMatrix,window,normopt);
     
    if(gtoSave == -1 && gfromSave == -1) {    // Si pas d'option -u et -v
      gto   = -1;
      gfrom = -1; 
      glen  = -1;
    }
    else {
      if(gfromSave != -1 && gtoSave != -1) {  // Si option -u && -v
	gto   = gtoSave;
	gfrom = gfromSave;
	glen  = -1;
      }
      else {
	if(gfromSave != -1) {                 // Si option -u
	  gto   = -1;
	  gfrom = gfromSave;
	  glen  = -1;
	}
	else {                                // Si option -v
	  gto   = gtoSave;
	  gfrom = -1;
	  glen  = -1;
	}
      }
    }
  }
  
  OpenDoor();
  
  if (!PorteOuverte && Data_Len > 6000) 
    {
      printf("Valid license key not found ! The software is therefore limited to 6kb long\n");
      printf("sequences. Please contact Thomas.Schiex@toulouse.inra.fr and ask for a FREE\n");
      printf("license key.\n");
      exit(2);
    }
  
  // ---------------------------------------------------------------------------
  /* Exploiting spliced alignements against EST and complete cDNA */
  // ---------------------------------------------------------------------------
  if (estopt) {
    FILE *fEST;
    strcpy(tempname,fstname);
    strcat(tempname,".est");
    fEST = FileOpen(NULL,tempname, "r");
    NumEST = 0;
    HitTable = ESTAnalyzer(fEST,ESTMatch,EstM,Don,Acc,&NumEST);
    fclose(fEST);
  }
  
  // ---------------------------------------------------------------------------
  // Lecture des RIKEN
  // ---------------------------------------------------------------------------
  if (raflopt) {
    FILE *fRAFL;
    int beg5,end5;
    int beg3,end3;
    char name[FILENAME_MAX+1];
    RAFLgene tmp;
    std::vector <RAFLgene> RAFLtmp;
    
    fprintf(stderr,"Reading RAFL gene...");
    fflush(stderr);
    strcpy(tempname,fstname);
    strcat(tempname,".riken");
    fRAFL = FileOpen(NULL,tempname, "r");
    while ((j=fscanf (fRAFL,"%d %d %*s %d %d %*s %s",
		      &beg5,&end5,&beg3,&end3,name)) == 5) {
      tmp.deb=Min(beg3,Min(end3,Min(beg5,end5)));
      tmp.fin=Max(beg3,Max(end3,Max(beg5,end5)));
      strcpy(tmp.ID,name);
      
      tmp.sens = (((beg5+end5) < (beg3+end3)) ?1:-1);
      // si le gene est trop court, on ne peut pas connaitre le sens !
      if (abs((beg5+end5) - (beg3+end3)) <100) tmp.sens = 0;
      
      RAFLtmp.push_back(tmp);
    }
    fclose(fRAFL);
    if (j != EOF) {
      fprintf(stderr,"Incorrect RAFL file\n");
      exit(2);
    }
    fprintf(stderr,"%d RAFL EST pairs read, ",RAFLtmp.size());
    
    sort(RAFLtmp.begin(),RAFLtmp.end(),Before);

    for (j =0; j<(int)RAFLtmp.size()-1 ;j++) {
      if (RAFLtmp[j].fin - RAFLtmp[j+1].deb >= 60) { // grand overlap
        fprintf(stderr,"fusion...");
        // grand overlap, on fusionne les deux riken (meme gene)
        RAFLtmp[j].deb = Min( RAFLtmp[j].deb , RAFLtmp[j+1].deb );
        RAFLtmp[j].fin = Max( RAFLtmp[j].fin , RAFLtmp[j+1].fin );
        if ( (RAFLtmp[j].sens==RAFLtmp[j+1].sens) || 
             ((RAFLtmp[j].sens==0) || (RAFLtmp[j+1].sens==0))) {
          RAFLtmp[j].sens= ( (RAFLtmp[j].sens==0) ? RAFLtmp[j+1].sens : RAFLtmp[j].sens );
          RAFL.push_back(RAFLtmp[j]); 
          // si pas de contradiction dans les sens, on prend le gene resultant 
        }
        j++;
        //on saute le suivant
      }
      else{
        if (RAFLtmp[j].fin - RAFLtmp[j+1].deb > 0) {
          fprintf(stderr,"overlap...");
          i= RAFLtmp[j].fin;
          RAFLtmp[j].fin= Min(i,RAFLtmp[j+1].deb);
          RAFLtmp[j+1].deb= Max(i,RAFLtmp[j+1].deb);
        }
        RAFL.push_back(RAFLtmp[j]);
      }
    }
    if(j==(int)RAFLtmp.size()-1) RAFL.push_back(RAFLtmp[j]);
    
    if (graph) {
      const int HLen = 30;
      
      for (j =0; j<(int)RAFL.size() ;j++) {
	PlotBarF(RAFL[j].deb, 0,0.9,0.2,2);
	PlotLine(RAFL[j].deb,RAFL[j].deb+HLen,0,0,1.0,1.0,2);
	PlotBarF(RAFL[j].fin, 0,0.9,0.2,2);
	PlotLine(RAFL[j].fin-HLen,RAFL[j].fin,0,0,1.0,1.0,2);
      }
    }
    
    fprintf(stderr,"%d kept\n",RAFL.size());
    fflush(stderr);
    if (RAFL.size() < 1) raflopt=FALSE;
  }
  
  // ---------------------------------------------------------------------------
  // Analysing NC input
  // ---------------------------------------------------------------------------
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
      for (i = deb; i <= fin; i++) {
	ForcedIG[i] = TRUE;
	if (graph) PlotBarI(i,0,0.25,2,6);
	}
    }
  }
    
  // ---------------------------------------------------------------------------
  // Blastx against protein databases 1st simple version, we simply
  // enhance coding probability according to phase Dangerous (NR
  // contains translation of frameshifted sequences) another
  // possibility would be too forbid introns/intergenic states
  // 8 levels of confidence may be used.
  // ---------------------------------------------------------------------------
  const int LevelColor[3] = {6,7,8}; 
  if (blastopt)
    {
      FILE *fblast;
      int  overlap,deb,fin,phase,ProtDeb,ProtFin,PProtFin,level;
      float score;
      int  Pfin = 0, PPhase = 0;
      char A[128],B[128];
      char *ProtId, *PProtId,*tmp;
      REAL GlobalScore;
      REAL PGlobalScore = 0;
      const int MaxOverlap = 10; 
      const int MaxHitLen  = 15000;
      
      fprintf(stderr,"Reading Blastx data, level... ");
      fflush(stderr);
      
      for( k=0; k<(int)strlen(blastArg); k++ )
	{
	  strcpy(tempname,fstname);
	  strcat(tempname,".blast");
	  i = strlen(tempname);
	  tempname[i] = blastArg[k] - 1;
	  tempname[i+1] = 0;
		 
	  fblast = FileOpen(NULL,tempname, "r");
	  if (fblast == NULL) continue;
	  
	  level = (blastArg[k] - '0') - 1;
	  A[0] = B[0]= 0;
	  ProtId = A;
	  PProtId = B;
	  PProtFin = -10000;
	  
	  while (fscanf(fblast,"%d %d %f %*s %d %s %d %d\n", 
			&deb, &fin, &score, &phase, ProtId,&ProtDeb,&ProtFin) != EOF) {
	    if (abs(fin-deb) > MaxHitLen) {
	      fprintf(stderr,"Similarity of extreme length rejected. Check %s\n",ProtId);
	      continue;
	    }
	    
	    if (phase < 0) {
	      j = deb;
	      deb = fin;
	      fin = j;
	      j = ProtDeb;
	      ProtDeb = ProtFin;
	      ProtFin = j;
	    }
	    GlobalScore=((REAL)score)/((REAL)abs(fin-deb));
	    
	    overlap=0;
	    // Reconstruction GAPS -> INTRONS
	    if ( (strcmp(ProtId,PProtId) == 0) &&
		 (abs(PProtFin-ProtDeb)<= (MaxOverlap)) ) {
	      // Detection d'un INTRON
	      overlap= (PProtFin+1-ProtDeb)*3; // *3 car coord.nucleique
	      // overlap >0 : hits chevauchants, overlap <0 : hits espaces
	      if ((deb-Pfin+overlap) >= MinLength[8] ){
		// Le tableau des introns est rempli prudemment: uniquement les bordures, et sans serrer pres de l'exon. 
		for (i= Pfin-(overlap<0)*overlap ; i< Pfin+MinLength[8]-abs(overlap) ; i++) {
		  // debut de l'intron seulement... (dont le score est fonction de l'exon precedent)
		  if (BlastS[level] >= ProtMatchLevel[i]){
		    if (BlastS[level] > ProtMatchLevel[i]){
		      ProtMatchLevel[i]= BlastS[level];
		      ProtMatch[i]= -PGlobalScore;
		      ProtMatchPhase[i]=0;
		    }
		    else{
		      if (PGlobalScore >= fabs(ProtMatch[i])){
			ProtMatch[i]= -PGlobalScore;
			ProtMatchPhase[i]= 0;
		      }
		    }
		  }
		  //		     j=((phase < 0)?-4:4);
		  //		     PlotBarI(i,j,0.6+(level/8.0),1,LevelColor[level]);
		}
		for (i = deb-MinLength[8]+abs(overlap) ; i < deb+(overlap<0)*overlap ; i++){
		  // ...et fin de l'intron (score est fonction de l'exon actuel)
		  if (BlastS[level] >= ProtMatchLevel[i]){
		    if (BlastS[level] > fabs(ProtMatchLevel[i])){
		      ProtMatchLevel[i]= BlastS[level];
		      ProtMatch[i]= -GlobalScore;
		      ProtMatchPhase[i]=0;
		    }
		    else {
		      if (GlobalScore >= fabs(ProtMatch[i])){
			ProtMatch[i]= -GlobalScore;
			ProtMatchPhase[i]= 0;
		      }
		    }
		  }
		  //		     j=((phase < 0)?-4:4);
		  //		     PlotBarI(i,j,0.6+(level/8.0),1,LevelColor[level]);
		}
	      }
	      if (graph && level<3) {
		PlotLine(Pfin,deb,PPhase,phase,0.6+(level/8.0),0.6+(level/8.0),LevelColor[level]);
		//		   for(i= Pfin-(overlap<0)*overlap; i < deb+(overlap<0)*overlap ; i++){
		//		     j=((phase < 0)?-4:4);
		//		     PlotBarI(i,j,0.6+(level/8.0),1,LevelColor[level]);
		//		   }
	      }
	    }
	    
	    PGlobalScore=GlobalScore;
	    Pfin = fin;
	    tmp = PProtId;
	    PProtId = ProtId;
	    ProtId = tmp;
	    PProtFin = ProtFin;
	    PPhase = phase;
	    
	    phase = ph06(phase);
	    
	    // HITS -> CODANT
	    if (graph && level<3) {
	      for (i = deb-1; i < fin; i++){
		PlotBarI(i,PhaseAdapt(phase),0.6+(level/8.0),1,LevelColor[level]);
	      }
	    }
	    
	    for (i = deb-1; i < fin; i++)  {
	      if (BlastS[level] >= ProtMatchLevel[i]){
		if (BlastS[level] > ProtMatchLevel[i]){
		  ProtMatchLevel[i] = BlastS[level];
		  ProtMatch[i]= GlobalScore;
		  ProtMatchPhase[i]= PhaseAdapt(phase);
		}
		else{
		  if (GlobalScore >= fabs(ProtMatch[i])){ 
		    ProtMatch[i] = GlobalScore;
		    ProtMatchPhase[i]= PhaseAdapt(phase);
		  }
		}
	      }
	    }
	  }
	  fprintf(stderr,"%d ",level);
	  fflush(stderr);
	  
	  fclose(fblast);       
	}
      fprintf(stderr," done\n");
    }
  // ---------------------------------------------------------------------------
  // Data allocation for the shortest path with length constraints
  // algorithm
  // ---------------------------------------------------------------------------
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
  // 17 Intergenique
  // ---------------------------------------------------------------------------
  
  char *Choice;
  BackPoint *LBP[18];
  REAL BestU;
  signed   char best   = 0;
  unsigned char Switch = 0;
  
  Choice =  new char[Data_Len+2];
  for (i = 0; i < Data_Len +2; i++) Choice[i] = 0;
  
  for (i = 0; i < 18; i++) {
    LBP[i] = new BackPoint(i, -1,0.0);
    LBP[i]->Next = LBP[i];
    LBP[i]->Prev = LBP[i];
  }
  
  // Les PrevBP sont des pointeurs sur les "opening edges"
  // Les PBest correspondent au cout du chemin correspondant
  
  REAL  maxi, PBest[26];
  BackPoint *PrevBP[26];
  int source = 0;

  // ---------------------------------------------------------------------------
  // Couts initiaux  
  // ---------------------------------------------------------------------------
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
  
  UserInfoList = ContentsUser;
  //  WriteUtils(UserInfoList,stdout);
  
  for (i = 0; i <= Data_Len; i++) {
    
    // compute coding... probabilities
    Fill_Score(TheSeq,IMMatrix,i,BaseScore);
    
    //application des regles utilisateur pour i
    // les Weight[] sont statiques ici (pas de signaux)
    if (userinfo) 
      Util(i,UserInfoList);
    
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
    
    // ---------------------------------------------------------------------------
    // Calcul de la position par rapport aux genes RAFL (Riken Ara.Full-Length)
    // valeurs de RAFLpos:
    // 0-> en dehors,
    // 1-> frontiere 5prime (intergenique ou UTR5  obligatoire),
    // 2-> frontiere 3prime (intergenique ou UTR3  obligatoire),
    // 3-> dedans(penalisation IG)
    // ---------------------------------------------------------------------------
    if (raflopt) {
      RAFLpos=0;
      if (i > RAFL[RAFLindex].fin) {
        // si on depasse le RAFL, on prend l'eventuel prochain
        if (RAFLindex+1 < (int)RAFL.size()) RAFLindex++; 
        else raflopt=FALSE;
      }

      // dans le Riken:
      if (i >= RAFL[RAFLindex].deb-2) {
        // bordure min:
        if (i == RAFL[RAFLindex].deb-2)
          RAFLpos= ( (RAFL[RAFLindex].sens== 1) ? 1 : 2);
        else {
          if (i < RAFL[RAFLindex].fin) RAFLpos=3;
          else {
            // frontiere max:
            if (i== RAFL[RAFLindex].fin) {
              RAFLpos= ( (RAFL[RAFLindex].sens==-1) ? 1 : 2);
            }
          }
        }
      }
    } 
    // ----------------------------------------------------------------
    // ------------------ Exons en forward ----------------------------
    // ----------------------------------------------------------------
    for (k = 0; k < 3; k++) {
      
      maxi = NINFINITY;     
      
      // S'il y  a un STOP en phase on ne peut continuer
      if ((i % 3 == k) && Stop[0][i]) 
	LBP[k]->Update(DontCrossStop); 
#ifdef PAYTOIGNORE      
      else // on ne prend pas le donneur
	LBP[k]->Update(log(1.0-Don[0][i]));
#endif
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
#ifndef PAYTOIGNORE
	  BestU -= log(1.0-Start[0][i]);
#endif 
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
#ifndef PAYTOIGNORE
      BestU -= log(1.0-Acc[0][i]);
#endif
      // Un test tordu pour casser le cou aux NaN
      if (isnan(maxi) || (BestU > maxi)) {
	maxi = BestU;
	best = 6+((i-k+3) % 3);
	Switch = SwitchAcc;
      }
      
      if (best != k) 
	LBP[k]->InsertNew(((best >= 18) ? source : best),Switch,i,maxi,PrevBP[best]);
      
      LBP[k]->Update(log(BaseScore[k])+log(4*TransitionWeight[0]));
       
      // Si on a un Gap EST ou si l'on connait le sens du match EST
      if (estopt)
	if ((ESTMatch[i] & Gap) || ((ESTMatch[i] & Hit) && !(ESTMatch[i] & HitForward)))
	  LBP[k]->Update(EstP);
      
      if (blastopt)
	if ((ProtMatch[i]<0) || ((ProtMatch[i]>0)&&(ProtMatchPhase[i]!=PhaseAdapt(k))))
	  LBP[k]->Update(-fabs(ProtMatch[i])*ProtMatchLevel[i]);
      
      if ((ForcedIG != NULL) && ForcedIG[i])
	LBP[k]->Update(IGPenalty);
      
      if (raflopt){
	if ((RAFLpos==1) || (RAFLpos==2) ||
	    ((RAFLpos==3) && (RAFL[RAFLindex].sens==-1)) )
	  LBP[k]->Update(RAFLPenalty);
      }
    }
    // ----------------------------------------------------------------
    // ------------------------- Exons en reverse ---------------------
    // ----------------------------------------------------------------
    for (k = 3; k<6; k++) {
      
      maxi = NINFINITY;
      
      // On continue sauf si l'on rencontre un autre STOP
      if (((Data_Len-i) % 3 == k-3) && Stop[1][i]) 
	LBP[k]->Update(DontCrossStop);
#ifdef PAYTOIGNORE      
      else // sinon on ne prend pas le donneur 
	LBP[k]->Update(log(1.0-Don[1][i]));
#endif
      
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
#ifndef PAYTOIGNORE
      BestU -= log(1.0-Don[1][i]);
#endif
      // Un test tordu pour casser le cou aux NaN
      if (isnan(maxi) || (BestU > maxi)) {
	maxi = BestU;
	best = 9+((Data_Len-i-k) % 3);
	Switch =  SwitchDon;
      }
            
      if (best != k) 
	LBP[k]->InsertNew(((best >= 19) ? source : best),Switch,i,maxi,PrevBP[best]);
      
      LBP[k]->Update(log(BaseScore[k])+log(4*TransitionWeight[0])) ;

      // Si on a un Gap EST ou si l'on connait le sens du match EST
      if (estopt)
	if ((ESTMatch[i] & Gap) || ((ESTMatch[i] & Hit) && !(ESTMatch[i] & HitReverse)))
	  LBP[k]->Update(EstP);

     if (blastopt)
       if ((ProtMatch[i]<0) || ((ProtMatch[i]>0)&&(ProtMatchPhase[i]!=PhaseAdapt(k))))
	LBP[k]->Update(-fabs(ProtMatch[i])*ProtMatchLevel[i]);

      if ((ForcedIG != NULL) && ForcedIG[i])
	LBP[k]->Update(IGPenalty);

      if (raflopt){
	if ((RAFLpos==1) || (RAFLpos==2) ||
	    ((RAFLpos==3) && (RAFL[RAFLindex].sens== 1)) )
	  LBP[k]->Update(RAFLPenalty);
      }
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
    
    LBP[InterGen5]->Update(log(BaseScore[8])+log(4*TransitionWeight[4]));
    //    LBP[InterGen5]->Update(log(TheSeq->GC_AT(i))+log(4));
    if (estopt)
      LBP[InterGen5]->Update(((ESTMatch[i] & (Gap|Hit)) != 0)*EstP);

    if(blastopt && (ProtMatch[i]!=0))
      LBP[InterGen5]->Update(-fabs(ProtMatch[i])*ProtMatchLevel[i]);

    if (raflopt){
      if (RAFLpos==3) LBP[InterGen5]->Update(RAFLPenalty);
    }

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
    
    LBP[InterGen3]->Update(log(BaseScore[8])+log(4*TransitionWeight[4]));
    //LBP[InterGen3]->Update(log(TheSeq->GC_AT(i))+log(4));
    if (estopt)
      LBP[InterGen3]->Update(((ESTMatch[i] & (Gap|Hit)) != 0)*EstP);

    if (blastopt && (ProtMatch[i] != 0))
      LBP[InterGen3]->Update(-fabs(ProtMatch[i])*ProtMatchLevel[i]);

    if (raflopt){
      if (RAFLpos==3) LBP[InterGen3]->Update(RAFLPenalty);
    }

    // ----------------------------------------------------------------
    // ---------------------- UTR 5' direct ---------------------------
    // ----------------------------------------------------------------
    maxi = NINFINITY;
    
#ifdef PAYTOIGNORE      
    // On reste 5' direct. On ne prend pas le Start eventuel.
    //  Kludge: si on a un EST qui nous dit que l'on est dans un
    //  intron, on oublie
    if (!estopt || (ESTMatch[i] & Gap) == 0)
      LBP[UTR5F]->Update(log(1.0-Start[0][i]));
#endif

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

    //    LBP[UTR5F]->Update(log(BaseScore[8])+log(3.999));
    //    LBP[UTR5F]->Update(log(TheSeq->GC_AT(i))+log(3.999));
    LBP[UTR5F]->Update(log(BaseScore[9])+log(4*TransitionWeight[2]));

    if (blastopt && (ProtMatch[i] != 0))
      LBP[UTR5F]->Update(-fabs(ProtMatch[i])*ProtMatchLevel[i]);

    if (raflopt){
      if ( (RAFLpos==2) ||
	   ((RAFL[RAFLindex].sens==-1) && (RAFLpos>=1) ) )
	LBP[UTR5F]->Update(RAFLPenalty);
    }
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

    //    LBP[UTR3F]->Update(log(BaseScore[9])+log(4));
    //    LBP[UTR3F]->Update(log(BaseScore[8])+log(4));
    LBP[UTR3F]->Update(log(BaseScore[11])+log(4*TransitionWeight[3]));
    //    LBP[UTR3F]->Update(log(TheSeq->GC_AT(i))+log(4));

    if (blastopt && (ProtMatch[i] != 0))
      LBP[UTR3F]->Update(-fabs(ProtMatch[i])*ProtMatchLevel[i]);

    if (raflopt){
      if ( (RAFLpos==1) ||
	   ((RAFL[RAFLindex].sens==-1) && (RAFLpos>=1) ) )
	LBP[UTR3F]->Update(RAFLPenalty);
    }

    // ----------------------------------------------------------------
    // ----------------------- UTR 5'reverse --------------------------
    // ----------------------------------------------------------------
    maxi = NINFINITY;
    
#ifdef PAYTOIGNORE      
    // On reste 5' reverse
    //  Kludge: si on a un EST qui nous dit que l'on est dans un
    //  intron, on oublie
    if (!estopt || (ESTMatch[i] & Gap) == 0)
      LBP[UTR5R]->Update(log(1.0-Start[1][i]));
#endif

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
#ifndef PAYTOIGNORE
	BestU -= log(1.0-Start[1][i]);
#endif
	// Un test tordu pour casser le cou aux NaN
	if (isnan(maxi) || (BestU > maxi)) {
	  maxi = BestU;
	  best = k+12;
	  Switch = SwitchStart;
	}
      }
    }

    if (best != -1) LBP[UTR5R]->InsertNew(best % 12,Switch,i,maxi,PrevBP[best]);

    //    LBP[UTR5R]->Update(log(BaseScore[8])+log(3.999));
    //    LBP[UTR5R]->Update(log(TheSeq->GC_AT(i))+log(3.999));
    LBP[UTR5R]->Update(log(BaseScore[10])+log(4*TransitionWeight[2]));

    if (blastopt && (ProtMatch[i]!=0))
      LBP[UTR5R]->Update(-fabs(ProtMatch[i])*ProtMatchLevel[i]);

    if (raflopt){
      if ( (RAFLpos==2) || 
	   ((RAFL[RAFLindex].sens== 1) && (RAFLpos>=1) ) )
	LBP[UTR5R]->Update(RAFLPenalty);
    }
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
      
    //    LBP[UTR3R]->Update(log(BaseScore[10])+log(4));
    //    LBP[UTR3R]->Update(log(BaseScore[8])+log(4));
    LBP[UTR3R]->Update(log(BaseScore[12])+log(4)*TransitionWeight[3]);
    //    LBP[UTR3R]->Update(log(TheSeq->GC_AT(i))+log(4));

    if (blastopt && (ProtMatch[i] != 0))
      LBP[UTR3R]->Update(-fabs(ProtMatch[i])*ProtMatchLevel[i]);

    if (raflopt){
      if ( (RAFLpos==1) ||
	   ((RAFL[RAFLindex].sens== 1) && (RAFLpos>=1)))
	LBP[UTR3R]->Update(RAFLPenalty);
    }

    // ----------------------------------------------------------------
    // ---------------- Introns de phase k forward --------------------
    // ----------------------------------------------------------------
    for (k = 0; k<3; k++) {

      maxi = NINFINITY;

      // On reste intronique
#ifdef PAYTOIGNORE      
      LBP[6+k]->Update(log(1.0-Acc[0][i]));
#endif
      LBP[6+k]->BestUsable(i,SwitchAny,0,&BestU);
      // Un test tordu pour casser le cou aux NaN
      if (isnan(maxi) || (BestU > maxi)) {
	maxi = BestU;
	best = -1;
      }
      // - on quitte un exon
      BestU = PBest[((i-k+3) % 3)]+log(Don[0][i]);
#ifndef PAYTOIGNORE
      BestU -= log(1.0-Don[0][i]);
#endif
      // Un test tordu pour casser le cou aux NaN
      if (isnan(maxi) || (BestU > maxi)) {
	maxi = BestU;
	best = ((i-k+3) % 3);
	Switch = SwitchDon;
      }
 
      if (best != -1) LBP[6+k]->InsertNew(best,Switch,i,maxi,PrevBP[best]);
      
      LBP[6+k]->Update(log(BaseScore[6])+log(4*TransitionWeight[1]));


      // Si on a un Hit EST ou si l'on connait le sens du match EST
      if (estopt)
	if((ESTMatch[i] & Hit) ||
	   ((ESTMatch[i] & Gap) && !(ESTMatch[i] & GapForward)))
	  LBP[6+k]->Update(EstP);

      if ((ForcedIG != NULL) && ForcedIG[i])
	LBP[6+k]->Update(IGPenalty);

      if (blastopt && (ProtMatch[i] > 0))
	LBP[6+k]->Update(-ProtMatch[i]*ProtMatchLevel[i]);

      if (raflopt){
	if ( (RAFLpos==1) || (RAFLpos==2) || 
		 ( (RAFLpos==3) && (RAFL[RAFLindex].sens==-1) ) )
	  LBP[6+k]->Update(RAFLPenalty);
      }
    }

    // ----------------------------------------------------------------
    // ----------------- Introns de phase -k reverse ------------------
    // ----------------------------------------------------------------
    for (k = 0; k<3; k++) {
      
      maxi = NINFINITY;
      
      // On reste intronique
#ifdef PAYTOIGNORE      
      LBP[9+k]->Update(log(1.0-Acc[1][i]));
#endif
      LBP[9+k]->BestUsable(i,SwitchAny,0,&BestU);
      // Un test tordu pour casser le cou aux NaN
      if (isnan(maxi) || (BestU > maxi)) {
	maxi = BestU;
	best = -1;
      }

      // On quitte un exon
      BestU = PBest[3+((Data_Len-i-k) % 3)]+log(Acc[1][i]);
#ifndef PAYTOIGNORE
      BestU -= log(1.0-Acc[1][i]);
#endif
      // Un test tordu pour casser le cou aux NaN
      if (isnan(maxi) || (BestU > maxi)) {
	maxi = BestU;
	best = 3+((Data_Len-i-k) % 3);
	Switch = SwitchAcc;
      }
      
      if (best != -1) LBP[9+k]->InsertNew(best,Switch,i,maxi,PrevBP[best]);
      
      LBP[9+k]->Update(log(BaseScore[7])+log(4*TransitionWeight[1]));

      // Si on a un Hit EST ou si l'on connait le sens du match EST
      if (estopt)
	if ((ESTMatch[i] & Hit) ||
      	  ((ESTMatch[i] & Gap) && !(ESTMatch[i] & GapReverse)))
      	LBP[9+k]->Update(EstP);

      if ((ForcedIG != NULL) && ForcedIG[i])
	LBP[9+k]->Update(IGPenalty);

      if (blastopt && (ProtMatch[i] > 0))
	LBP[9+k]->Update(-ProtMatch[i]*ProtMatchLevel[i]);

      if (raflopt){
	if ( (RAFLpos==1) || (RAFLpos==2) || 
		 ( (RAFLpos==3) && (RAFL[RAFLindex].sens== 1)) )
	  LBP[9+k]->Update(RAFLPenalty);
      }
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
  LBP[UTR5F]->Update(log(FivePrimePrior/2.0)/2.0);
  LBP[UTR3F]->Update(log(ThreePrimePrior/2.0)/2.0);  
  LBP[UTR5R]->Update(log(FivePrimePrior/2.0)/2.0);
  LBP[UTR3R]->Update(log(ThreePrimePrior/2.0)/2.0);

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
  
  fprintf(stderr,"Optimal path length = %#f\n",-maxi+log(4)*(Data_Len+1));

  // Sanity check ! A feasible path has not been found ?
  if (isnan(maxi))
    fprintf(stderr,"WARNING: no feasible path, inconsistent data !\n");

  if (graph) PlotPredictions(Data_Len,Choice,Stop,Start,Acc,Don);
  
#include "Output.h"
  // free used memory
  if (graph) {
    fprintf(stderr,"Dumping images (\"%s.---.png\")...",grname);
    fflush(stderr);
    ClosePNG();
    fprintf(stderr, "done\n\n");
  }
  delete TheSeq;
  
  if (HitTable) {
    delete HitTable[NumEST];
    delete HitTable;
  }
  
  HitTable = NULL;
  
  delete [] Stop[0];
  delete [] Stop[1];

  if (estopt)  delete [] ESTMatch;

  if (blastopt) {
    delete [] ProtMatch;
    delete [] ProtMatchLevel;
    delete [] ProtMatchPhase;
  }

  delete [] Start[0];
  delete [] Start[1];

  delete [] Don[0];
  delete [] Don[1];
  
  delete [] Acc[0];
  delete [] Acc[1];

  delete [] Choice;
  } // fin de traitement de chaque séquence....

  // free remaining used memory
  if (IMMatrix[5] != IMMatrix[3]) delete  IMMatrix[5];
  if (IMMatrix[6] != IMMatrix[3]) delete  IMMatrix[6];

  for  (i = 0;  i < 5;  i ++)
    delete  IMMatrix[i];

  return  0;
}
