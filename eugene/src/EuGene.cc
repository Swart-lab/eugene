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
// remettre Frameshifts
// alignement cds cDNA et proteines sur les splice sites

// supprimer l'approximation longueur single (on peut rater un bon
// epissage si un meilleur START l'occulte mais qu'il est trop pres
// pour etre utilisable). Il faudrait faire 3 pistes single + 3 pistes
// non single.

#define PAYTOIGNORE 1

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

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
#include <string.h>

#ifdef HAVE_STRINGS_H
#include <strings.h>
#endif

//#include "dmalloc.h"
//#include "dmallocc.cc"
//#define STAND 1

#include "Const.h"
#include "SensorIF.h"
#include "MSensor.h"
#include "Param.h"
#include "DNASeq.h"
#include "BackP.h"
#include "clef.h"
#include "Output.h"
#include "Prediction.h"

#ifndef FILENAME_MAX
#define FILENAME_MAX        1024
#endif

// ------------------ Globals --------------------
MasterSensor MS;
Parameters   PAR;

// -------------------------------------------------------------------------
//            MAIN
// -------------------------------------------------------------------------
int main  (int argc, char * argv [])
{
  int    i, j, k;
  FILE   *fp;
  char   *grnameFile = NULL;
  char   grname[FILENAME_MAX+1];
  DATA	 Data;
  DNASeq *TheSeq;
  int    Data_Len;
  double ExPrior, InPrior, IGPrior, FivePrior, ThreePrior;
  int    minL[18], minDiv, minFlow, minConv, min5, min3;
  int    graph;
  const unsigned char SwitchMask[18] = 
    {SwitchAny,SwitchAny,SwitchAny,SwitchAny,SwitchAny,SwitchAny,
     SwitchAny,SwitchAny,SwitchAny,SwitchAny,SwitchAny,SwitchAny,
     SwitchStart,SwitchStart,SwitchStart,SwitchStop,SwitchStop,SwitchStop};
 
  // Lecture de la ligne d'arg et du fichier .par
  PAR.initParam(argc, argv);
  
  // Objectif -> limiter le nombre d'appel à la map de PAR
  ExPrior    = PAR.getD("EuGene.ExonPrior");
  InPrior    = PAR.getD("EuGene.IntronPrior");
  IGPrior    = PAR.getD("EuGene.InterPrior"); 
  FivePrior  = PAR.getD("EuGene.FivePrimePrior");
  ThreePrior = PAR.getD("EuGene.ThreePrimePrior");
  minL[0]  = PAR.getI("EuGene.minL0");
  minL[1]  = PAR.getI("EuGene.minL1");
  minL[2]  = PAR.getI("EuGene.minL2");
  minL[3]  = PAR.getI("EuGene.minL3");
  minL[4]  = PAR.getI("EuGene.minL4");
  minL[5]  = PAR.getI("EuGene.minL5");
  minL[6]  = PAR.getI("EuGene.minL6");
  minL[7]  = PAR.getI("EuGene.minL7");
  minL[8]  = PAR.getI("EuGene.minL8");
  minL[9]  = PAR.getI("EuGene.minL9");
  minL[10] = PAR.getI("EuGene.minL10");
  minL[11] = PAR.getI("EuGene.minL11");
  minL[12] = PAR.getI("EuGene.minL12");
  minL[13] = PAR.getI("EuGene.minL13");
  minL[14] = PAR.getI("EuGene.minL14");
  minL[15] = PAR.getI("EuGene.minL15");
  minL[16] = PAR.getI("EuGene.minL16");
  minL[17] = PAR.getI("EuGene.minL17");
  minDiv  = PAR.getI("EuGene.minDiv");
  minFlow = PAR.getI("EuGene.minFlow");
  minConv = PAR.getI("EuGene.minConv");
  min5    = PAR.getI("EuGene.min5Prime");
  min3    = PAR.getI("EuGene.min3Prime");
  graph = PAR.getI("Output.graph");
  
  MS.InitMaster();
  
  ReadKey(PAR.getC("EuGene.key"), "EUGENEAT");
  
  // ---------------------------------------------------------------------------  
  // Lecture de la sequence
  // ---------------------------------------------------------------------------  
  int sequence;
  for (sequence = optind; sequence < argc ; sequence++) {
    
    PAR.set("fstname", argv[sequence]);

    // read fasta file
    fp = (*PAR.getC("fstname") ? FileOpen (NULL,PAR.getC("fstname"), "r") : stdin);
    
    if (fp == NULL) {
      fprintf(stderr, "cannot open fasta file %s\n", PAR.getC("fstname"));
      exit(3);
    }
    
    fprintf(stderr,"-----------------------------------");
    fprintf(stderr,"--------------------------------\nLoading sequence...");
    fflush(stderr);
    
    TheSeq = new DNASeq(fp);
    
    if (fp != stdin) fclose(fp);
    
    Data_Len = TheSeq->SeqLen;
    
    fprintf(stderr,"%s, %d bases read, ",TheSeq->Name, Data_Len);
    
    fprintf (stderr,"GC Proportion = %.1f%%\n", 
	     (TheSeq->Markov0[BitG]+TheSeq->Markov0[BitC])*100.0);
    
        
    // ---------------------------------------------------------------------------
    // Preparation sortie graphique + Scores
    // ---------------------------------------------------------------------------
    if (graph) {
      int gto       = PAR.getI("Output.gto");
      int gfrom     = PAR.getI("Output.gfrom");
      int glen      = PAR.getI("Output.glen");

      // Récupération du nom du fichier d'origine
      grnameFile = BaseName(PAR.getC("fstname"));
      
      // Construction du nom de sortie (*.png)
      strcpy(grname, PAR.getC("Output.Prefix"));
      strcat(grname, grnameFile);
      *rindex(grname, '.') = 0;             // on enleve l'extension (.fasta typ.)
      if(PAR.getC("grnameArg")[0] != '0') { // -p h sans -g ou -g NO pas d'argument
	strcat(grname,".");
        strcat(grname,PAR.getC("grnameArg"));
      }

      if ((gfrom <= 0)|| (gfrom >= Data_Len))
	gfrom = 1;
      if ((gto <= 0)  || (gto <= gfrom) || (gto > Data_Len))
	gto = Data_Len;
      
      if ((PAR.getI("Output.gfrom") != -1) || PAR.getI("Output.gto") != -1) {
	sprintf(grname+strlen(grname), ".%d", gfrom);
	sprintf(grname+strlen(grname), "-%d", gto);
      }
      
      gfrom--;
      gto--;
      
      if (glen < 0)
	glen = ((gto-gfrom+1 <= 6000) ? gto-gfrom+1 : 6000);
      
      InitPNG(PAR.getI("Output.resx"), PAR.getI("Output.resy"), PAR.getI("Output.offset"),
	      gfrom, gto, PAR.getI("Output.golap"), glen, grname);
    }
    
    OpenDoor();
    
    if (!PorteOuverte && Data_Len > 6000) 
      {
	printf("Valid license key not found ! The software is therefore limited to 6kb long\n");
	printf("sequences. Please contact Thomas.Schiex@toulouse.inra.fr and ask for a FREE\n");
	printf("license key.\n");
	exit(2);
      }
    
    MS.InitSensors(TheSeq);

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
    Prediction *pred;
    BackPoint  *LBP[18];
    REAL BestU;
    double NormalizingPath = 0.0;
    signed   char best   = 0;
    unsigned char Switch = 0;
    
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
    LBP[ExonF1]->Update(log(ExPrior/6.0)/2.0);
    LBP[ExonF2]->Update(log(ExPrior/6.0)/2.0);
    LBP[ExonF3]->Update(log(ExPrior/6.0)/2.0);
    LBP[ExonR1]->Update(log(ExPrior/6.0)/2.0);
    LBP[ExonR2]->Update(log(ExPrior/6.0)/2.0);
    LBP[ExonR3]->Update(log(ExPrior/6.0)/2.0);
    
    // Intron
    LBP[IntronF1]->Update(log(InPrior/6.0)/2.0);
    LBP[IntronF2]->Update(log(InPrior/6.0)/2.0);
    LBP[IntronF3]->Update(log(InPrior/6.0)/2.0);
    LBP[IntronR1]->Update(log(InPrior/6.0)/2.0);
    LBP[IntronR2]->Update(log(InPrior/6.0)/2.0);
    LBP[IntronR3]->Update(log(InPrior/6.0)/2.0);
    
    // Intergenique 
    LBP[InterGen5]->Update(log(IGPrior)/2.0); // ameliorations sur 5' reverse
    LBP[InterGen3]->Update(log(IGPrior)/2.0); // ameliorations sur 3' direct
    
    // UTR 5' et 3'
    LBP[UTR5F]->Update(log(FivePrior /2.0)/2.0);
    LBP[UTR3F]->Update(log(ThreePrior/2.0)/2.0);  
    LBP[UTR5R]->Update(log(FivePrior /2.0)/2.0);
    LBP[UTR3R]->Update(log(ThreePrior/2.0)/2.0);
    
    for (i = 0; i <= Data_Len; i++) {

    MS.GetInfoAt(TheSeq, i, &Data);
      
      maxi = -NINFINITY;
      for (k = 0 ; k < 18; k++) {
	LBP[k]->BestUsable(i,SwitchAny,0,&BestU);
	if ((BestU > NINFINITY) && (BestU < maxi)) maxi =BestU;
     }

      for (k = 0 ; k < 18; k++) 
	LBP[k]->Update(-maxi);
      NormalizingPath += maxi;

      // Calcul des meilleures opening edges
      for(k=0; k<18; k++)
	PrevBP[k] = LBP[k%12]->BestUsable(i, SwitchMask[k], minL[k], &PBest[k]);
      
      // intergenic: longueur min depend du sens
      // -> -> ou <- <-   MinFlow
      // -> <-            MinConv
      // <- ->            MinDiv
      PrevBP[23] = LBP[InterGen5]->StrictBestUsable(i, minDiv,  &PBest[23]);
      PrevBP[24] = LBP[InterGen5]->StrictBestUsable(i, minFlow, &PBest[24]);
      
      PrevBP[18] = LBP[InterGen3]->StrictBestUsable(i, minFlow, &PBest[18]);
      PrevBP[25] = LBP[InterGen3]->StrictBestUsable(i, minConv, &PBest[25]);
      
      // UTR 5' et 3' direct
      PrevBP[19] = LBP[UTR5F]->StrictBestUsable(i, min5, &PBest[19]);
      PrevBP[20] = LBP[UTR3F]->StrictBestUsable(i, min3, &PBest[20]);
      
      // UTR 5' et 3' reverse
      PrevBP[21] = LBP[UTR5R]->StrictBestUsable(i, min5, &PBest[21]);
      PrevBP[22] = LBP[UTR3R]->StrictBestUsable(i, min3, &PBest[22]);
      
      // ----------------------------------------------------------------
      // ------------------ Exons en forward ----------------------------
      // ----------------------------------------------------------------
      for (k = 0; k < 3; k++) {
	maxi = NINFINITY;     
	
	// On va tout droit.
	// S'il y  a un STOP en phase on ne peut continuer
	if ((i % 3 == k))
	  LBP[k]->Update(Data.sig[DATA::Stop].weight[Signal::ForwardNo]);
#ifdef PAYTOIGNORE      
	else // on ne prend pas le donneur
	  LBP[k]->Update(Data.sig[DATA::Don].weight[Signal::ForwardNo]);
#endif
	LBP[k]->BestUsable(i,SwitchAny,0,&BestU);
	// Un test tordu pour casser le cou aux NaN
	if (isnan(maxi) || (BestU > maxi)) {
	  maxi = BestU;
	  best = k;
	}
	
	// On commence a coder (Start)
	// Ca vient d'une UTR 5' forward
	if ((i % 3 == k)) {
	  BestU = PBest[19]+Data.sig[DATA::Start].weight[Signal::Forward];
#ifndef PAYTOIGNORE
	  BestU -= Data.sig[DATA::Start].weight[Signal::ForwardNo];
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
	BestU = PBest[6+((i-k+3) % 3)]+Data.sig[DATA::Acc].weight[Signal::Forward];
#ifndef PAYTOIGNORE
	BestU -= Data.sig[DATA::Acc].weight[Signal::ForwardNo];
#endif
	// Un test tordu pour casser le cou aux NaN
	if (isnan(maxi) || (BestU > maxi)) {
	  maxi = BestU;
	  best = 6+((i-k+3) % 3);
	  Switch = SwitchAcc;
	}
	
	// Il y a une insertion (frameshift). Pour l'instant, on ne
	// prend pas en compte le saut de nucléotide.
	BestU = PBest[(k+2)%3] + Data.sig[DATA::Ins].weight[Signal::Forward];
	if (isnan(maxi) || (BestU > maxi)) {
	  maxi = BestU;
	  best = (k+2)%3;
	  Switch = SwitchIns;
	}

	// Il y a une deletion (frameshift)
	BestU = PBest[(k+1)%3] + Data.sig[DATA::Del].weight[Signal::Forward];
	if (isnan(maxi) || (BestU > maxi)) {
	  maxi = BestU;
	  best = (k+1)%3;
	  Switch = SwitchDel;
	}

	if (best != k) 
	  LBP[k]->InsertNew(((best >= 18) ? source : best),Switch,i,maxi,PrevBP[best]);
	LBP[k]->Update(Data.contents[k]);
      }
      
      // ----------------------------------------------------------------
      // ------------------------- Exons en reverse ---------------------
      // ----------------------------------------------------------------
      for (k = 3; k<6; k++) {
	maxi = NINFINITY;
	
	// On continue sauf si l'on rencontre un autre STOP
	if (((Data_Len-i) % 3 == k-3)) 
	  LBP[k]->Update(Data.sig[DATA::Stop].weight[Signal::ReverseNo]);
#ifdef PAYTOIGNORE      
	else // sinon on ne prend pas le donneur 
	  LBP[k]->Update(Data.sig[DATA::Don].weight[Signal::ReverseNo]);
#endif
	
	LBP[k]->BestUsable(i,SwitchAny,0,&BestU);
	// Un test tordu pour casser le cou aux NaN
	if (isnan(maxi) || (BestU > maxi)) {
	  maxi = BestU;
	  best = k;
	}
	
	// On commence a coder (Stop)
	// Ca vient d'une UTR 3' reverse
	if (((Data_Len-i) % 3 == k-3)) {
	  BestU = PBest[22] + Data.sig[DATA::Stop].weight[Signal::Reverse];
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
	BestU = PBest[9+((Data_Len-i-k) % 3)] + Data.sig[DATA::Don].weight[Signal::Reverse];
#ifndef PAYTOIGNORE
	BestU -= Data.sig[DATA::Don].weight[Signal::ReverseNo];
#endif
	// Un test tordu pour casser le cou aux NaN
	if (isnan(maxi) || (BestU > maxi)) {
	  maxi = BestU;
	  best = 9+((Data_Len-i-k) % 3);
	  Switch =  SwitchDon;
	}

	// Il y a une insertion (frameshift)
	BestU = PBest[3+(k+2)%3] + Data.sig[DATA::Ins].weight[Signal::Reverse];
	if (isnan(maxi) || (BestU > maxi)) {
	  maxi = BestU;
	  best = 3+(k+2)%3;
	  Switch = SwitchIns;
	}

	// Il y a une deletion (frameshift)
	BestU = PBest[3+(k+1)%3] + Data.sig[DATA::Del].weight[Signal::Reverse];
	if (isnan(maxi) || (BestU > maxi)) {
	  maxi = BestU;
	  best = 3+(k+1)%3;
	  Switch = SwitchDel;
	}

	if (best != k) 
	  LBP[k]->InsertNew(((best >= 19) ? source : best),Switch,i,maxi,PrevBP[best]);
	LBP[k]->Update(Data.contents[k]);
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
      BestU = PBest[21] + Data.sig[DATA::tStart].weight[Signal::Reverse];
      // Un test tordu pour casser le cou aux NaN
      if (isnan(maxi) || (BestU > maxi)) {
	maxi = BestU;
	best = 21;
	source  = 15;
	Switch = SwitchTransStart;
      }
      
      if (best != -1) 
	LBP[InterGen5]->InsertNew(source,Switch,i,maxi,PrevBP[best]);
      
      LBP[InterGen5]->Update(Data.contents[8]);
      
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
      BestU = PBest[20] + Data.sig[DATA::tStop].weight[Signal::Forward];
      // Un test tordu pour casser le cou aux NaN
      if (isnan(maxi) || (BestU > maxi)) {
	maxi = BestU;
	best = 20;
	source  = 14;
	Switch = SwitchTransStop;
      }
      
      if (best != -1) 
	LBP[InterGen3]->InsertNew(source,Switch,i,maxi,PrevBP[best]);
      
      LBP[InterGen3]->Update(Data.contents[8]);
      
      // ----------------------------------------------------------------
      // ---------------------- UTR 5' direct ---------------------------
      // ----------------------------------------------------------------
      maxi = NINFINITY;
      
#ifdef PAYTOIGNORE
      // On reste 5' direct. On ne prend pas le Start eventuel.
      //  Kludge: si on a un EST qui nous dit que l'on est dans un
      //  intron, on oublie
      if (!PAR.getI("Sensor.Est.use") || (Data.ESTMATCH_TMP & Gap) == 0)  // WARNING
	LBP[UTR5F]->Update(Data.sig[DATA::Start].weight[Signal::ForwardNo]);
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
      BestU = PBest[23] + Data.sig[DATA::tStart].weight[Signal::Forward];
      // Un test tordu pour casser le cou aux NaN
      if (isnan(maxi) || (BestU > maxi)) {
	maxi = BestU;
	best = 23;
	source = InterGen5;
	Switch = SwitchTransStart;
      }
      
      // Sur 3' direct
      BestU = PBest[18] + Data.sig[DATA::tStart].weight[Signal::Forward];
      // Un test tordu pour casser le cou aux NaN
      if (isnan(maxi) || (BestU > maxi)) {
	maxi = BestU;
	best = 18;
	source = InterGen3;
	Switch = SwitchTransStart;
      }
      
      if (best != -1) LBP[UTR5F]->InsertNew(source,Switch,i,maxi,PrevBP[best]);

      LBP[UTR5F]->Update(Data.contents[9]);
      
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
	if ((i % 3 == k)) {
	  BestU = PBest[k+12] + Data.sig[DATA::Stop].weight[Signal::Forward];
	  // Un test tordu pour casser le cou aux NaN
	  if (isnan(maxi) || (BestU > maxi)) {
	    maxi = BestU;
	    best = k+12;
	    Switch = SwitchStop;
	  }
	}
      }
      
      if (best != -1) LBP[UTR3F]->InsertNew(best %12,Switch,i,maxi,PrevBP[best]);
      
      LBP[UTR3F]->Update(Data.contents[11]);
      
      // ----------------------------------------------------------------
      // ----------------------- UTR 5'reverse --------------------------
      // ----------------------------------------------------------------
      maxi = NINFINITY;
      
#ifdef PAYTOIGNORE      
      // On reste 5' reverse
      //  Kludge: si on a un EST qui nous dit que l'on est dans un
      //  intron, on oublie
      if (!PAR.getI("Sensor.Est.use") || (Data.ESTMATCH_TMP & Gap) == 0)  // WARNING
	LBP[UTR5R]->Update(Data.sig[DATA::Start].weight[Signal::ReverseNo]);
#endif
      
      LBP[UTR5R]->BestUsable(i,SwitchAny,0,&BestU);
      // Un test tordu pour casser le cou aux NaN
      if (isnan(maxi) || (BestU > maxi)) {
	maxi = BestU;
	best = -1;
      }
      
      // Ca vient d'un exon reverse + START
      for (k = 3; k < 6; k++) {
	if (((Data_Len-i) % 3 == k-3)) {
	  BestU = PBest[k+12]+Data.sig[DATA::Start].weight[Signal::Reverse];
#ifndef PAYTOIGNORE
	  BestU -= Data.sig[DATA::Start].weight[Signal::ReverseNo];
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
      
      LBP[UTR5R]->Update(Data.contents[10]);
      
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
      BestU = PBest[24] + Data.sig[DATA::tStop].weight[Signal::Reverse];
      // Un test tordu pour casser le cou aux NaN
      if (isnan(maxi) || (BestU > maxi)) {
	maxi = BestU;
	best = 24;
	source = 12;
	Switch = SwitchTransStop;
      }
      
      // Sur 3' direct
      BestU = PBest[25] + Data.sig[DATA::tStop].weight[Signal::Reverse];
      // Un test tordu pour casser le cou aux NaN
      if (isnan(maxi) || (BestU > maxi)) {
	maxi = BestU;
	best = 25;
	source = 12;
	Switch = SwitchTransStop;
      }
      
      if (best != -1) LBP[UTR3R]->InsertNew(source,Switch,i,maxi,PrevBP[best]);

      LBP[UTR3R]->Update(Data.contents[12]);
      
      // ----------------------------------------------------------------
      // ---------------- Introns de phase k forward --------------------
      // ----------------------------------------------------------------
      for (k = 0; k<3; k++) {
	maxi = NINFINITY;
	
	// On reste intronique
#ifdef PAYTOIGNORE      
	LBP[6+k]->Update(Data.sig[DATA::Acc].weight[Signal::ForwardNo]);
#endif
	LBP[6+k]->BestUsable(i,SwitchAny,0,&BestU);
	// Un test tordu pour casser le cou aux NaN
	if (isnan(maxi) || (BestU > maxi)) {
	  maxi = BestU;
	  best = -1;
	}
	// - on quitte un exon
	BestU = PBest[((i-k+3) % 3)]+ Data.sig[DATA::Don].weight[Signal::Forward];
#ifndef PAYTOIGNORE
	BestU -= Data.sig[DATA::Don].weight[Signal::ForwardNo];
#endif
	// Un test tordu pour casser le cou aux NaN
	if (isnan(maxi) || (BestU > maxi)) {
	  maxi = BestU;
	  best = ((i-k+3) % 3);
	  Switch = SwitchDon;
	}
	
	if (best != -1) LBP[6+k]->InsertNew(best,Switch,i,maxi,PrevBP[best]);

	LBP[6+k]->Update(Data.contents[6]);
      }
      
      // ----------------------------------------------------------------
      // ----------------- Introns de phase -k reverse ------------------
      // ----------------------------------------------------------------
      for (k = 0; k<3; k++) {
	maxi = NINFINITY;
	
	// On reste intronique
#ifdef PAYTOIGNORE      
	LBP[9+k]->Update(Data.sig[DATA::Acc].weight[Signal::ReverseNo]);
#endif
	LBP[9+k]->BestUsable(i,SwitchAny,0,&BestU);
	// Un test tordu pour casser le cou aux NaN
	if (isnan(maxi) || (BestU > maxi)) {
	  maxi = BestU;
	  best = -1;
	}
	
	// On quitte un exon
	BestU = PBest[3+((Data_Len-i-k) % 3)]+Data.sig[DATA::Acc].weight[Signal::Reverse];
#ifndef PAYTOIGNORE
	BestU -= Data.sig[DATA::Acc].weight[Signal::ReverseNo];
#endif
	// Un test tordu pour casser le cou aux NaN
	if (isnan(maxi) || (BestU > maxi)) {
	  maxi = BestU;
	  best = 3+((Data_Len-i-k) % 3);
	  Switch = SwitchAcc;
	}
	
	if (best != -1) LBP[9+k]->InsertNew(best,Switch,i,maxi,PrevBP[best]);
	
	LBP[9+k]->Update(Data.contents[7]);
      }
    }
    
    LBP[ExonF1]->Update(log(ExPrior/6.0)/2.0);
    LBP[ExonF2]->Update(log(ExPrior/6.0)/2.0);
    LBP[ExonF3]->Update(log(ExPrior/6.0)/2.0);
    LBP[ExonR1]->Update(log(ExPrior/6.0)/2.0);
    LBP[ExonR2]->Update(log(ExPrior/6.0)/2.0);
    LBP[ExonR3]->Update(log(ExPrior/6.0)/2.0);

    LBP[IntronF1]->Update(log(InPrior/6.0)/2.0);
    LBP[IntronF2]->Update(log(InPrior/6.0)/2.0);
    LBP[IntronF3]->Update(log(InPrior/6.0)/2.0);
    LBP[IntronR1]->Update(log(InPrior/6.0)/2.0);
    LBP[IntronR2]->Update(log(InPrior/6.0)/2.0);
    LBP[IntronR3]->Update(log(InPrior/6.0)/2.0);
    
    // Intergenique 
    LBP[InterGen5]->Update(log(IGPrior)/2.0); 
    LBP[InterGen3]->Update(log(IGPrior)/2.0); 
    
    // UTR 5' et 3'
    LBP[UTR5F]->Update(log(FivePrior /2.0)/2.0);
    LBP[UTR3F]->Update(log(ThreePrior/2.0)/2.0);  
    LBP[UTR5R]->Update(log(FivePrior /2.0)/2.0);
    LBP[UTR3R]->Update(log(ThreePrior/2.0)/2.0);
    
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
        
    pred = LBP[j]->BackTrace();
        
    for  (i = 0;  i < 18;  i ++) LBP[i]->Zap();
    
    if (!PorteOuverte && Data_Len > 6000) 
      exit(2);
    
    fprintf(stderr,"Optimal path length = %.4f\n",-maxi-NormalizingPath);
    
    // Sanity check ! A feasible path has not been found ?
    if (isnan(maxi+NormalizingPath) || maxi+NormalizingPath == NINFINITY)
      fprintf(stderr,"WARNING: no feasible path, inconsistent data !\n");
    
    if (graph)
      pred->plotPred();
    
    Output(TheSeq, pred, sequence, argc, argv);
    MS.PostAnalyse(pred);
    
    // Free used memory
    if (graph) {
      fprintf(stderr,"\nDumping images (\"%s.---.png\")...", grname);
      fflush(stderr);
      ClosePNG();
      fprintf(stderr, "done\n");
    }

    delete TheSeq;
    pred->resetPred();
    //    MS.ResetType();

  }// fin de traitement de chaque séquence....
 
  MS.ResetSensors();

  fprintf(stderr,"-----------------------------------");
  fprintf(stderr,"--------------------------------\n");
  
  return  0;
}
