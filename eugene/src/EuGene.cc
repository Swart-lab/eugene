//-----------------------------------------------------------------------------
//  This program finds exons/introns and intergenic regions (including UTR)
//  Copyright T. Schiex 1999
//
//  $Id$
//-----------------------------------------------------------------------------
// Manage T/TG/TA for introns in coding regions.

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

#include "Const.h"
#include "SensorIF.h"
#include "MSensor.h"
#include "Param.h"
#include "DNASeq.h"
#include "BackP.h"
#include "PenaltyDist.h"
#include "Output.h"
#include "Prediction.h"
#include "../Parametrization/ParaOptimization.h"

#ifndef FILENAME_MAX
#define FILENAME_MAX        1024
#endif

// ------------------ Globals --------------------
MasterSensor*    MS;
Parameters       PAR;
ParaOptimization OPTIM;

// -------------------------------------------------------------------------
// Predict a gene
// -------------------------------------------------------------------------
Prediction* Predict (DNASeq* TheSeq, MasterSensor* MSensor)
{

  int    j, k;
  double ExPrior, InPrior, IGPrior, FivePrior, ThreePrior,IntronFivePrior;
  double PredictTime = cpuTime();
  int    Data_Len;
  DATA	 Data;
  
  Prediction *pred;
  Track  LBP[NbTracks];
  double BestU;
  double NormalizingPath = 0.0;
  signed   char best   = 0;
  
  // Objectif -> limiter le nombre d'appel à la map de PAR
  ExPrior    = PAR.getD("EuGene.ExonPrior");
  InPrior    = PAR.getD("EuGene.IntronPrior");
  IGPrior    = PAR.getD("EuGene.InterPrior"); 
  FivePrior  = PAR.getD("EuGene.FivePrimePrior");
  ThreePrior = PAR.getD("EuGene.ThreePrimePrior");
  IntronFivePrior = InPrior;

  Data_Len = TheSeq->SeqLen;
  
  for (j = 0; j < NbTracks; j++) {
    LBP[j].Path.SetState(j);
  }

  LBP[InitF1].LoadPenalty(PAR.getC("EuGene.InitExDist"));
  LBP[InitF2].LoadPenalty(PAR.getC("EuGene.InitExDist"));
  LBP[InitF3].LoadPenalty(PAR.getC("EuGene.InitExDist"));
  LBP[InitR1].LoadPenalty(PAR.getC("EuGene.InitExDist"));
  LBP[InitR2].LoadPenalty(PAR.getC("EuGene.InitExDist"));
  LBP[InitR3].LoadPenalty(PAR.getC("EuGene.InitExDist"));

  LBP[IntronF1].LoadPenalty(PAR.getC("EuGene.IntronDist"));
  LBP[IntronF2].LoadPenalty(PAR.getC("EuGene.IntronDist"));
  LBP[IntronF3].LoadPenalty(PAR.getC("EuGene.IntronDist"));
  LBP[IntronR1].LoadPenalty(PAR.getC("EuGene.IntronDist"));
  LBP[IntronR2].LoadPenalty(PAR.getC("EuGene.IntronDist"));
  LBP[IntronR3].LoadPenalty(PAR.getC("EuGene.IntronDist"));
  
  LBP[SnglF1].LoadPenalty(PAR.getC("EuGene.SnglExDist"));
  LBP[SnglF2].LoadPenalty(PAR.getC("EuGene.SnglExDist"));
  LBP[SnglF3].LoadPenalty(PAR.getC("EuGene.SnglExDist"));
  LBP[SnglR1].LoadPenalty(PAR.getC("EuGene.SnglExDist"));
  LBP[SnglR2].LoadPenalty(PAR.getC("EuGene.SnglExDist"));
  LBP[SnglR3].LoadPenalty(PAR.getC("EuGene.SnglExDist"));

  LBP[IntrF1].LoadPenalty(PAR.getC("EuGene.IntrExDist"));
  LBP[IntrF2].LoadPenalty(PAR.getC("EuGene.IntrExDist"));
  LBP[IntrF3].LoadPenalty(PAR.getC("EuGene.IntrExDist"));
  LBP[IntrR1].LoadPenalty(PAR.getC("EuGene.IntrExDist"));
  LBP[IntrR2].LoadPenalty(PAR.getC("EuGene.IntrExDist"));
  LBP[IntrR3].LoadPenalty(PAR.getC("EuGene.IntrExDist"));

  LBP[TermF1].LoadPenalty(PAR.getC("EuGene.TermExDist"));
  LBP[TermF2].LoadPenalty(PAR.getC("EuGene.TermExDist"));
  LBP[TermF3].LoadPenalty(PAR.getC("EuGene.TermExDist"));
  LBP[TermR1].LoadPenalty(PAR.getC("EuGene.TermExDist"));
  LBP[TermR2].LoadPenalty(PAR.getC("EuGene.TermExDist"));
  LBP[TermR3].LoadPenalty(PAR.getC("EuGene.TermExDist"));

  LBP[InterGen].LoadPenalty(PAR.getC("EuGene.InterGDist"));
  
  // UTR 5' et 3'
  LBP[UTR5F].LoadPenalty(PAR.getC("EuGene.5PrimeDist"));
  LBP[UTR3F].LoadPenalty(PAR.getC("EuGene.3PrimeDist"));
  LBP[UTR5R].LoadPenalty(PAR.getC("EuGene.5PrimeDist"));
  LBP[UTR3R].LoadPenalty(PAR.getC("EuGene.3PrimeDist"));

  // Introns d'UTR5
  LBP[IntronU5F].LoadPenalty(PAR.getC("EuGene.IntronDist"));
  LBP[IntronU5R].LoadPenalty(PAR.getC("EuGene.IntronDist"));

  // Introns d'UTR3
  LBP[IntronU3F].LoadPenalty(PAR.getC("EuGene.IntronDist"));
  LBP[IntronU3R].LoadPenalty(PAR.getC("EuGene.IntronDist"));

  // Les PrevBP sont des pointeurs sur les "opening edges"
  // Les PBest correspondent au cout du chemin correspondant
  double  maxi, PBest[NbTracks];
  BackPoint *PrevBP[NbTracks];
  
  // --------------------------------------------------------------------------
  // Couts initiaux  
  // --------------------------------------------------------------------------
  // Codant
  LBP[InitF1].Update(log(ExPrior/6.0)/2.0);
  LBP[InitF2].Update(log(ExPrior/6.0)/2.0);
  LBP[InitF3].Update(log(ExPrior/6.0)/2.0);
  LBP[InitR1].Update(log(ExPrior/6.0)/2.0);
  LBP[InitR2].Update(log(ExPrior/6.0)/2.0);
  LBP[InitR3].Update(log(ExPrior/6.0)/2.0);
  
  // Intron
  LBP[IntronF1].Update(log(InPrior/6.0)/2.0);
  LBP[IntronF2].Update(log(InPrior/6.0)/2.0);
  LBP[IntronF3].Update(log(InPrior/6.0)/2.0);
  LBP[IntronR1].Update(log(InPrior/6.0)/2.0);
  LBP[IntronR2].Update(log(InPrior/6.0)/2.0);
  LBP[IntronR3].Update(log(InPrior/6.0)/2.0);
  
  // Sngl 
  LBP[SnglF1].Update(log(ExPrior/6.0)/2.0);
  LBP[SnglF2].Update(log(ExPrior/6.0)/2.0);
  LBP[SnglF3].Update(log(ExPrior/6.0)/2.0);
  LBP[SnglR1].Update(log(ExPrior/6.0)/2.0);
  LBP[SnglR2].Update(log(ExPrior/6.0)/2.0);
  LBP[SnglR3].Update(log(ExPrior/6.0)/2.0);

  // Intr 
  LBP[IntrF1].Update(log(ExPrior/6.0)/2.0);
  LBP[IntrF2].Update(log(ExPrior/6.0)/2.0);
  LBP[IntrF3].Update(log(ExPrior/6.0)/2.0);
  LBP[IntrR1].Update(log(ExPrior/6.0)/2.0);
  LBP[IntrR2].Update(log(ExPrior/6.0)/2.0);
  LBP[IntrR3].Update(log(ExPrior/6.0)/2.0);

  // Term
  LBP[TermF1].Update(log(ExPrior/6.0)/2.0);
  LBP[TermF2].Update(log(ExPrior/6.0)/2.0);
  LBP[TermF3].Update(log(ExPrior/6.0)/2.0);
  LBP[TermR1].Update(log(ExPrior/6.0)/2.0);
  LBP[TermR2].Update(log(ExPrior/6.0)/2.0);
  LBP[TermR3].Update(log(ExPrior/6.0)/2.0);

  // Intergenique 
  LBP[InterGen].Update(log(IGPrior)/2.0); 
  
  // UTR 5' et 3'
  LBP[UTR5F].Update(log(FivePrior /2.0)/2.0);
  LBP[UTR3F].Update(log(ThreePrior/2.0)/2.0);  
  LBP[UTR5R].Update(log(FivePrior /2.0)/2.0);
  LBP[UTR3R].Update(log(ThreePrior/2.0)/2.0);

  // Introns d'UTR5
  LBP[IntronU5F].Update(log(IntronFivePrior/2.0)/2.0);
  LBP[IntronU5R].Update(log(IntronFivePrior/2.0)/2.0);

  // Introns d'UTR5
  LBP[IntronU3F].Update(log(IntronFivePrior/2.0)/2.0);
  LBP[IntronU3R].Update(log(IntronFivePrior/2.0)/2.0);

  // --------------------------------------------------------------------------
  // Demarrage de la programmation dynamique
  // --------------------------------------------------------------------------
  for (int nuc = 0; nuc <= Data_Len; nuc++) {

    // recuperation des infos
    MSensor->GetInfoAt(TheSeq, nuc, &Data);
    
    // normalisation des couts 
    maxi = -NINFINITY;
    for (k = 0 ; k < NbTracks; k++) {
      BestU = LBP[k].Optimal;
      if ((BestU > NINFINITY) && (BestU < maxi)) maxi =BestU;
    }

    for (k = 0 ; k < NbTracks; k++) 
      LBP[k].Update(-maxi);
    NormalizingPath += maxi;

#define ISPOSSIBLE(X,Y) (!(isinf(Data.sig[DATA::X].weight[Signal::Y])))
#define INEED(K) if (!PrevBP[K]) PrevBP[K] = LBP[K].BestUsable(nuc, &PBest[K])
#define PICOMP(C,S,B,O) if (C && ISPOSSIBLE(S,B)) {\
	BestU = PBest[O]+Data.sig[DATA::S].weight[Signal::B];\
	if (isnan(maxi) || (BestU > maxi)) {maxi = BestU; best = O;}}
#define INSERTANDPAYTHESLOPE(P,C)					\
    LBP[P].InsertNew(best, nuc, maxi,PrevBP[best]);\
    if (nuc < Data_Len){\
      LBP[P].Update(Data.contents[C]);\
      LBP[P].PayTheSlope();\
    }

    // ----------------------------------------------------------------
    // Calcul des meilleures opening edges
    // ----------------------------------------------------------------
    for (k = 0; k < NbTracks; k++) {
      PrevBP[k] = NULL;
      PBest[k] = NINFINITY;
    }
    
    // Exons F (de tous poils)
    for (k = 0; k <= 3; k++) {
      if ((nuc % 3 == k) && ISPOSSIBLE(Start,Forward))             INEED(UTR5F);
      if (ISPOSSIBLE(Acc,Forward))                                 INEED(IntronF1+((nuc-k+3)%3));
      if (ISPOSSIBLE(Ins,Forward)) { 	                           INEED(InitF1+(k+2)%3);
								   INEED(IntrF1+(k+2)%3);
      								   INEED(TermF1+(k+2)%3); }
      if (ISPOSSIBLE(Del,Forward)) {	                           INEED(InitF1+(k+1)%3);
      								   INEED(IntrF1+(k+1)%3);
      								   INEED(TermF1+(k+1)%3);
      }
    }

    // Exons R  (de tout poil)
    for (k = 3; k <= 6; k++) {
      if (((Data_Len-nuc) % 3 == k-3) && ISPOSSIBLE(Stop,Reverse)) INEED(UTR3R);
      if (ISPOSSIBLE(Don,Reverse))				   INEED(IntronR1+((Data_Len-nuc-k) % 3));
      if (ISPOSSIBLE(Ins,Reverse)) {				   INEED(InitR1+(k+2)%3);
								   INEED(IntrR1+(k+2)%3);
							  	   INEED(TermR1+(k+2)%3); }
      if (ISPOSSIBLE(Del,Reverse)) {				   INEED(InitR1+(k+1)%3);
								   INEED(IntrR1+(k+1)%3);
								   INEED(TermR1+(k+1)%3); }
    }

    // Introns F et R
    for (k = 0; k <= 3; k++) {
      if (ISPOSSIBLE(Don,Forward)) {				   INEED(InitF1+((nuc-k+3) % 3));
								   INEED(IntrF1+((nuc-k+3) % 3)); }
      if (ISPOSSIBLE(Acc,Reverse)) {				   INEED(IntrR1+((Data_Len-nuc-k) % 3));
      								   INEED(TermR1+((Data_Len-nuc-k) % 3)); }
    }
    
    // Intergenique
    if (ISPOSSIBLE(tStart,Reverse))     			   INEED(UTR5R);
    if (ISPOSSIBLE(tStop,Forward))			           INEED(UTR3F);

    // UTR 5' direct
    if (ISPOSSIBLE(tStart,Forward))			           INEED(InterGen);
    if (ISPOSSIBLE(Acc,Forward))			           INEED(IntronU5F);

    // Introns d'UTR5F
    if (ISPOSSIBLE(Don,Forward))                                   INEED(UTR5F);

    // Introns d'UTR3F
    if (ISPOSSIBLE(Don,Forward))                                   INEED(UTR3F);

    // UTR 3' direct
    if (ISPOSSIBLE(Acc,Forward))			           INEED(IntronU3F);
    if (ISPOSSIBLE(Stop,Forward)) {			           INEED(TermF1+nuc%3);
			   			                   INEED(SnglF1+nuc%3); }
    // UTR 5'reverse
    if (ISPOSSIBLE(Start,Reverse)) {                               INEED(InitR1+((Data_Len-nuc) % 3));
    				                                   INEED(SnglR1+((Data_Len-nuc) % 3)); }
    if (ISPOSSIBLE(Don,Reverse))                                   INEED(IntronU5R);

    // Introns d'UTR5R
    if (ISPOSSIBLE(Acc,Reverse))                                   INEED(UTR5R);

    // Introns d'UTR5R
    if (ISPOSSIBLE(Acc,Reverse))                                   INEED(UTR3R);

    // UTR 3' reverse
    if (ISPOSSIBLE(tStop,Reverse))                                 INEED(InterGen);
    if (ISPOSSIBLE(Don,Reverse))                                   INEED(IntronU3R);

    // ----------------------------------------------------------------
    // ------------------ Inits en forward ----------------------------
    // ----------------------------------------------------------------
    for (k = 0; k < 3; k++) {
      maxi = NINFINITY;     
      
      // On commence a coder (Start),si en phase ca vient d'une UTR 5' forward
      PICOMP((nuc % 3 == k),Start,Forward,UTR5F);
      // Il y a une insertion (frameshift). Saut de nucléotide ignore.
      PICOMP(true,Ins,Forward,InitF1+(k+2)%3);
      // Il y a une deletion (frameshift)
      PICOMP(true,Del,Forward,InitF1+(k+1)%3);
    
      // On va tout droit.
      // S'il y  a un STOP en phase on ne peut continuer
      if ((nuc % 3 == k))
	LBP[InitF1+k].Update(Data.sig[DATA::Stop].weight[Signal::ForwardNo]);
      LBP[InitF1+k].Update(Data.sig[DATA::Don].weight[Signal::ForwardNo]);

      INSERTANDPAYTHESLOPE(InitF1+k,k);
    }
    // ----------------------------------------------------------------
    // ------------------------- Inits en reverse ---------------------
    // ----------------------------------------------------------------
    for (k = 3; k<6; k++) {
      maxi = NINFINITY;
      
      // - on recommence a coder (Donneur) Ca vient d'un intron
      PICOMP(true,Don,Reverse,IntronR1+((Data_Len-nuc-k) % 3));
      // Il y a une insertion (frameshift)
      PICOMP(true,Ins,Reverse, InitR1+(k+2)%3);
      // Il y a une deletion (frameshift)
      PICOMP(true,Del,Reverse, InitR1+(k+1)%3);

      if (((Data_Len-nuc) % 3 == k-3)) 
	LBP[InitF1+k].Update(Data.sig[DATA::Stop].weight[Signal::ReverseNo]);
      LBP[InitF1+k].Update(Data.sig[DATA::Don].weight[Signal::ReverseNo]);

      INSERTANDPAYTHESLOPE(InitF1+k,k);
    }
    // ----------------------------------------------------------------
    // ------------------ Sngl en forward ----------------------------
    // ----------------------------------------------------------------
    for (k = 0; k < 3; k++) {
      maxi = NINFINITY;     
      
      // On commence a coder (Start),si en phase ca vient d'une UTR 5' forward
      PICOMP((nuc % 3 == k),Start,Forward,UTR5F);
      // Il y a une insertion (frameshift). Saut de nucléotide ignore.
      PICOMP(true,Ins,Forward,SnglF1+(k+2)%3);
      // Il y a une deletion (frameshift)
      PICOMP(true,Del,Forward,SnglF1+(k+1)%3);
    
      // On va tout droit.
      // S'il y  a un STOP en phase on ne peut continuer
      if ((nuc % 3 == k))
	LBP[SnglF1+k].Update(Data.sig[DATA::Stop].weight[Signal::ForwardNo]);
      LBP[SnglF1+k].Update(Data.sig[DATA::Don].weight[Signal::ForwardNo]);
      
      INSERTANDPAYTHESLOPE(SnglF1+k,k);
    }
    // ----------------------------------------------------------------
    // ------------------------- Sngl en reverse ---------------------
    // ----------------------------------------------------------------
    for (k = 3; k<6; k++) {
      maxi = NINFINITY;
      
      // On commence a coder (Stop). Ca vient d'une UTR 3' reverse
      PICOMP(((Data_Len-nuc) % 3 == k-3),Stop,Reverse,UTR3R);
      // Il y a une insertion (frameshift)
      PICOMP(true,Ins,Reverse, SnglR1+(k+2)%3);
      // Il y a une deletion (frameshift)
      PICOMP(true,Del,Reverse, SnglR1+(k+1)%3);

      if (((Data_Len-nuc) % 3 == k-3)) 
	LBP[SnglF1+k].Update(Data.sig[DATA::Stop].weight[Signal::ReverseNo]);
      LBP[SnglF1+k].Update(Data.sig[DATA::Don].weight[Signal::ReverseNo]);

      INSERTANDPAYTHESLOPE(SnglF1+k,k);
    }
    // ----------------------------------------------------------------
    // ------------------ Intrs en forward ----------------------------
    // ----------------------------------------------------------------
    for (k = 0; k < 3; k++) {
      maxi = NINFINITY;     
      
      // On recommence a coder (Accepteur). Ca vient d'un intron
      PICOMP(true,Acc,Forward,IntronF1+((nuc-k+3) % 3));
      // Il y a une insertion (frameshift). Saut de nucléotide ignore.
      PICOMP(true,Ins,Forward,IntrF1+(k+2)%3);
      // Il y a une deletion (frameshift)
      PICOMP(true,Del,Forward,IntrF1+(k+1)%3);
    
      // On va tout droit.
      // S'il y  a un STOP en phase on ne peut continuer
      if ((nuc % 3 == k))
	LBP[IntrF1+k].Update(Data.sig[DATA::Stop].weight[Signal::ForwardNo]);
      LBP[IntrF1+k].Update(Data.sig[DATA::Don].weight[Signal::ForwardNo]);

      INSERTANDPAYTHESLOPE(IntrF1+k,k);
    }
    // ----------------------------------------------------------------
    // ------------------------- Intrs en reverse ---------------------
    // ----------------------------------------------------------------
    for (k = 3; k<6; k++) {
      maxi = NINFINITY;
      
      // - on recommence a coder (Donneur) Ca vient d'un intron
      PICOMP(true,Don,Reverse,IntronR1+((Data_Len-nuc-k) % 3));
      // Il y a une insertion (frameshift)
      PICOMP(true,Ins,Reverse, IntrR1+(k+2)%3);
      // Il y a une deletion (frameshift)
      PICOMP(true,Del,Reverse, IntrR1+(k+1)%3);

      if (((Data_Len-nuc) % 3 == k-3)) 
	LBP[IntrF1+k].Update(Data.sig[DATA::Stop].weight[Signal::ReverseNo]);
      LBP[IntrF1+k].Update(Data.sig[DATA::Don].weight[Signal::ReverseNo]);

      INSERTANDPAYTHESLOPE(IntrF1+k,k);
    }
    // ----------------------------------------------------------------
    // ------------------ Terms en forward ----------------------------
    // ----------------------------------------------------------------
    for (k = 0; k < 3; k++) {
      maxi = NINFINITY;     
      
      // On recommence a coder (Accepteur). Ca vient d'un intron
      PICOMP(true,Acc,Forward,IntronF1+((nuc-k+3) % 3));
      // Il y a une insertion (frameshift). Saut de nucléotide ignore.
      PICOMP(true,Ins,Forward,TermF1+(k+2)%3);
      // Il y a une deletion (frameshift)
      PICOMP(true,Del,Forward,TermF1+(k+1)%3);
    
      // On va tout droit.
      // S'il y  a un STOP en phase on ne peut continuer
      if ((nuc % 3 == k))
	LBP[TermF1+k].Update(Data.sig[DATA::Stop].weight[Signal::ForwardNo]);
      LBP[TermF1+k].Update(Data.sig[DATA::Don].weight[Signal::ForwardNo]);

      INSERTANDPAYTHESLOPE(TermF1+k,k);
    }
    // ----------------------------------------------------------------
    // ------------------------- Terms en reverse ---------------------
    // ----------------------------------------------------------------
    for (k = 3; k<6; k++) {
      maxi = NINFINITY;
      
      // On commence a coder (Stop). Ca vient d'une UTR 3' reverse
      PICOMP(((Data_Len-nuc) % 3 == k-3),Stop,Reverse,UTR3R);
      // Il y a une insertion (frameshift)
      PICOMP(true,Ins,Reverse, TermR1+(k+2)%3);
      // Il y a une deletion (frameshift)
      PICOMP(true,Del,Reverse, TermR1+(k+1)%3);

      if (((Data_Len-nuc) % 3 == k-3)) 
	LBP[TermF1+k].Update(Data.sig[DATA::Stop].weight[Signal::ReverseNo]);
      LBP[TermF1+k].Update(Data.sig[DATA::Don].weight[Signal::ReverseNo]);

      INSERTANDPAYTHESLOPE(TermF1+k,k);
    }
    // ----------------------------------------------------------------
    // ------------------------ Intergenique --------------------------
    // ----------------------------------------------------------------
    // Ca peut venir d'une fin de 3' direct ou de 5' reverse
    maxi = NINFINITY;
    
    // From 5' reverse
    PICOMP(true,tStart,Reverse, UTR5R);
    // From 3' direct
    PICOMP(true,tStop,Forward, UTR3F);
 
    // On reste intergenique
    // et les transstartNO/TransstopNO ???
    LBP[InterGen].InsertNew(best, nuc, maxi,PrevBP[best]);

    INSERTANDPAYTHESLOPE(InterGen,DATA::InterG);
    // ----------------------------------------------------------------
    // ---------------------- UTR 5' direct ---------------------------
    // ----------------------------------------------------------------
    maxi = NINFINITY;
    
    // On vient de l'intergenique. 
    PICOMP(true,tStart,Forward, InterGen);
    // On peut venir aussi d'un intron d'UTR.
    PICOMP(true,Acc,Forward, IntronU5F);

    // On reste 5' direct. On ne prend pas le Start eventuel.
    LBP[UTR5F].Update(Data.sig[DATA::Start].weight[Signal::ForwardNo]);

    INSERTANDPAYTHESLOPE(UTR5F,DATA::UTR5F);
    // ----------------------------------------------------------------
    // ------------------- Introns d'UTR5F ----------------------------
    // ----------------------------------------------------------------
    maxi = NINFINITY;
      
    // on quitte l'UTR 
    PICOMP(true,Don,Forward, UTR5F);

    // On reste intronique
    LBP[IntronU5F].Update(Data.sig[DATA::Acc].weight[Signal::ForwardNo]);

    INSERTANDPAYTHESLOPE(IntronU5F,DATA::IntronUTRF);
    // ----------------------------------------------------------------
    // ------------------- Introns d'UTR3F ----------------------------
    // ----------------------------------------------------------------
    maxi = NINFINITY;
      
    // on quitte l'UTR 
    PICOMP(true,Don,Forward, UTR3F);

    // On reste intronique
    LBP[IntronU3F].Update(Data.sig[DATA::Acc].weight[Signal::ForwardNo]);

    INSERTANDPAYTHESLOPE(IntronU3F,DATA::IntronUTRF);
    // ----------------------------------------------------------------
    // ---------------------- UTR 3' direct ---------------------------
    // ----------------------------------------------------------------
    maxi = NINFINITY;

    // Ca vient d'un Term ou d'un Sngl direct + STOP
    PICOMP(true,Stop,Forward, TermF1+nuc%3);
    PICOMP(true,Stop,Forward, SnglF1+nuc%3);
    // On peut venir aussi d'un intron d'UTR.
    PICOMP(true,Acc,Forward, IntronU3F);

    // On reste 3' direct
    LBP[UTR3F].InsertNew(best, nuc, maxi,PrevBP[best]);

    INSERTANDPAYTHESLOPE(UTR3F,DATA::UTR3F);
    // ----------------------------------------------------------------
    // ----------------------- UTR 5'reverse --------------------------
    // ----------------------------------------------------------------
    maxi = NINFINITY;
    
    // Ca vient d'un Sngl ou Init reverse + START
    PICOMP(true,Start,Reverse, InitR1+((Data_Len-nuc) % 3));
    PICOMP(true,Start,Reverse, SnglR1+((Data_Len-nuc) % 3));
    // On peut venir aussi d'un intron d'UTR.
    PICOMP(true,Don,Reverse, IntronU5R);

    // On reste 5' reverse
    LBP[UTR5R].Update(Data.sig[DATA::Start].weight[Signal::ReverseNo]);

    INSERTANDPAYTHESLOPE(UTR5R,DATA::UTR5R);
    // ----------------------------------------------------------------
    // ------------------- Introns d'UTR5R ----------------------------
    // ----------------------------------------------------------------
    maxi = NINFINITY;
    
    // On quitte une UTR5R
    PICOMP(true,Acc,Reverse, UTR5R);

    // On reste intronique
    LBP[IntronU5R].Update(Data.sig[DATA::Acc].weight[Signal::ReverseNo]);

    INSERTANDPAYTHESLOPE(IntronU5R,DATA::IntronUTRR);
    // ----------------------------------------------------------------
    // ------------------- Introns d'UTR3R ----------------------------
    // ----------------------------------------------------------------
    maxi = NINFINITY;
    
    // On quitte une UTR5R
    PICOMP(true,Acc,Reverse, UTR3R);

    // On reste intronique
    LBP[IntronU3R].Update(Data.sig[DATA::Acc].weight[Signal::ReverseNo]);

    INSERTANDPAYTHESLOPE(IntronU3R,DATA::IntronUTRR);
    // ----------------------------------------------------------------
    // ----------------------- UTR 3'reverse --------------------------
    // ----------------------------------------------------------------
    maxi = NINFINITY;
    
    // On demarre depuis l'intergenique
    PICOMP(true,tStop,Reverse, InterGen);
    // On peut venir aussi d'un intron d'UTR.
    PICOMP(true,Don,Reverse, IntronU3R);

    // On reste 3' reverse
    LBP[UTR3R].InsertNew(best, nuc, maxi,PrevBP[best]);

    INSERTANDPAYTHESLOPE(UTR3R,DATA::UTR3R);
    // ----------------------------------------------------------------
    // ---------------- Introns de phase k forward --------------------
    // ----------------------------------------------------------------
    for (k = 0; k<3; k++) {
      maxi = NINFINITY;
      
      // - on quitte un Init ou un Intr
      PICOMP(true,Don,Forward, InitF1+((nuc-k+3) % 3));
      PICOMP(true,Don,Forward, IntrF1+((nuc-k+3) % 3));
      
      // On reste intronique
      LBP[IntronF1+k].Update(Data.sig[DATA::Acc].weight[Signal::ForwardNo]);
      
      INSERTANDPAYTHESLOPE(IntronF1+k,DATA::IntronF);
    }
    // ----------------------------------------------------------------
    // ----------------- Introns de phase -k reverse ------------------
    // ----------------------------------------------------------------
    for (k = 0; k<3; k++) {
      maxi = NINFINITY;
      
      // On quitte un Intr ou un Term reverse
      PICOMP(true,Acc,Reverse, IntrR1+((Data_Len-nuc-k) % 3));
      PICOMP(true,Acc,Reverse, TermR1+((Data_Len-nuc-k) % 3));

      // On reste intronique
      LBP[IntronR1+k].Update(Data.sig[DATA::Acc].weight[Signal::ReverseNo]);

      INSERTANDPAYTHESLOPE(IntronR1+k,DATA::IntronR);
    }
  }
  
  LBP[InitF1].Update(log(ExPrior/6.0)/2.0);
  LBP[InitF2].Update(log(ExPrior/6.0)/2.0);
  LBP[InitF3].Update(log(ExPrior/6.0)/2.0);
  LBP[InitR1].Update(log(ExPrior/6.0)/2.0);
  LBP[InitR2].Update(log(ExPrior/6.0)/2.0);
  LBP[InitR3].Update(log(ExPrior/6.0)/2.0);
  
  LBP[IntronF1].Update(log(InPrior/6.0)/2.0);
  LBP[IntronF2].Update(log(InPrior/6.0)/2.0);
  LBP[IntronF3].Update(log(InPrior/6.0)/2.0);
  LBP[IntronR1].Update(log(InPrior/6.0)/2.0);
  LBP[IntronR2].Update(log(InPrior/6.0)/2.0);
  LBP[IntronR3].Update(log(InPrior/6.0)/2.0);

  LBP[SnglF1].Update(log(ExPrior/6.0)/2.0);
  LBP[SnglF2].Update(log(ExPrior/6.0)/2.0);
  LBP[SnglF3].Update(log(ExPrior/6.0)/2.0);
  LBP[SnglR1].Update(log(ExPrior/6.0)/2.0);
  LBP[SnglR2].Update(log(ExPrior/6.0)/2.0);
  LBP[SnglR3].Update(log(ExPrior/6.0)/2.0);
 
  LBP[IntrF1].Update(log(ExPrior/6.0)/2.0);
  LBP[IntrF2].Update(log(ExPrior/6.0)/2.0);
  LBP[IntrF3].Update(log(ExPrior/6.0)/2.0);
  LBP[IntrR1].Update(log(ExPrior/6.0)/2.0);
  LBP[IntrR2].Update(log(ExPrior/6.0)/2.0);
  LBP[IntrR3].Update(log(ExPrior/6.0)/2.0);
 
  LBP[TermF1].Update(log(ExPrior/6.0)/2.0);
  LBP[TermF2].Update(log(ExPrior/6.0)/2.0);
  LBP[TermF3].Update(log(ExPrior/6.0)/2.0);
  LBP[TermR1].Update(log(ExPrior/6.0)/2.0);
  LBP[TermR2].Update(log(ExPrior/6.0)/2.0);
  LBP[TermR3].Update(log(ExPrior/6.0)/2.0);
 
  // Intergenique 
  LBP[InterGen].Update(log(IGPrior)/2.0); 
  
  // UTR 5' et 3'
  LBP[UTR5F].Update(log(FivePrior /2.0)/2.0);
  LBP[UTR3F].Update(log(ThreePrior/2.0)/2.0);  
  LBP[UTR5R].Update(log(FivePrior /2.0)/2.0);
  LBP[UTR3R].Update(log(ThreePrior/2.0)/2.0);

  // Introns d'UTR5
  LBP[IntronU5F].Update(log(IntronFivePrior/2.0)/2.0);
  LBP[IntronU5R].Update(log(IntronFivePrior/2.0)/2.0);
  // Introns d'UTR3
  LBP[IntronU3F].Update(log(IntronFivePrior/2.0)/2.0);
  LBP[IntronU3R].Update(log(IntronFivePrior/2.0)/2.0);
  
  for (j = 0; j < NbTracks; j++) {
    PrevBP[j] = LBP[j].BestUsable(Data_Len+1,&PBest[j],0);
    LBP[j].ForceNew(j,Data_Len+1,PBest[j],PrevBP[j]);
  }
  
  // Dump the Paths
//   for  (i = 0;  i < NbTracks;  i ++) {
//     printf("Piste %d\n",i);
//     LBP[i].Dump();
//   }

  // Backtrace in the graph
  maxi = PBest[0];
  j = 0;
  
  for (k = 1; k < NbTracks ; k++) {
    BestU = PBest[k];
    // Un test tordu pour casser le cou aux NaN
    if (isnan(maxi) || (BestU > maxi)) {
      maxi = BestU;
      j = k;
    }
  }

  pred = LBP[j].BackTrace();
  
  pred->OptimalPath = maxi+NormalizingPath;
  // Sanity check ! A feasible path has not been found ?
  if (isnan(pred->OptimalPath) || pred->OptimalPath == NINFINITY) {
    fprintf(stderr,"WARNING: no feasible path, inconsistent data !\n");
    pred->resetPred();
  }

  // Memory usage
  j=0;
  for (k =0;k<NbTracks;k++)  j+= LBP[k].NumBP;
  fprintf(stderr,"Number of Backpoints allocated: %d (%.1f sec)\n",
	  j,cpuTime()-PredictTime);

  // clean memory
  for  (k = 0;  k < NbTracks;  k++) LBP[k].Zap();

  return pred;
}

// -------------------------------------------------------------------------
// Read a fasta file
// -------------------------------------------------------------------------
DNASeq* ReadSequence (char* sequence_name)
{
  DNASeq *TheSeq;
  
  fprintf(stderr,"-------------------------------------");
  fprintf(stderr,"--------------------------------\nLoading sequence...");
  fflush(stderr);
  
  TheSeq = new DNASeq(sequence_name);
  
  fprintf(stderr,"%s, %d bases read, ",TheSeq->Name, TheSeq->SeqLen);
  
  fprintf (stderr,"GC Proportion = %.1f%%\n", 
	   (TheSeq->Markov0[BitG]+TheSeq->Markov0[BitC])*100.0);
  
  return TheSeq;
}

// -------------------------------------------------------------------------
//            MAIN
// -------------------------------------------------------------------------
int main  (int argc, char * argv [])
{
  DNASeq *TheSeq;
  int    Data_Len;
  Prediction *pred;
  char   *grnameFile = NULL;
  char   grname[FILENAME_MAX+1];
  int    graph;

  fprintf(stderr,"-------------------------------------");
  fprintf(stderr,"--------------------------------\n");

  // Lecture de la ligne d'arg et du fichier .par
  PAR.initParam(argc, argv);

  if ((std::string) PAR.getC("ParaOptimization.Use") ==  "TRUE") 
    OPTIM.ParaOptimize(argc, argv);
  else 
    {
      // Objectif : limiter les appels à la MAP
      graph = PAR.getI("Output.graph");
      
      int sequence;
      for (sequence = optind; sequence < argc ; sequence++) {
	
	PAR.set("fstname", argv[sequence]);
	// --------------------------------------------------------------------
	// Lecture de la sequence    
	// --------------------------------------------------------------------
	TheSeq = ReadSequence( PAR.getC("fstname") );
	Data_Len = TheSeq->SeqLen;
    
	// --------------------------------------------------------------------
	// Init MasterSensor
	// --------------------------------------------------------------------
	MS = new MasterSensor();
	MS->InitMaster(TheSeq);

	// --------------------------------------------------------------------
	// Preparation sortie graphique + Scores
	// --------------------------------------------------------------------
	if (graph) {
	  int gto       = PAR.getI("Output.gto");
	  int gfrom     = PAR.getI("Output.gfrom");
	  int glen      = PAR.getI("Output.glen");
      
	  // Récupération du nom du fichier d'origine
	  grnameFile = BaseName(PAR.getC("fstname"));
	  
	  // Construction du nom de sortie (*.png)
	  strcpy(grname, PAR.getC("Output.Prefix"));
	  strcat(grname, grnameFile);
	  *rindex(grname, '.') = 0;           // on enleve l'extension
	  if(PAR.getC("grnameArg")[0]!='0') { // -ph sans -g ou -g NO pas d'arg
	    strcat(grname,".");
	    strcat(grname,PAR.getC("grnameArg"));
	  }
	  
	  if ((gfrom <= 0)|| (gfrom >= Data_Len))
	    gfrom = 1;
	  if ((gto <= 0)  || (gto <= gfrom) || (gto > Data_Len))
	    gto = Data_Len;
	  
	  if ((PAR.getI("Output.gfrom")!=-1) || PAR.getI("Output.gto")!=-1) {
	    sprintf(grname+strlen(grname), ".%d", gfrom);
	    sprintf(grname+strlen(grname), "-%d", gto);
	  }
	  
	  gfrom--;
	  gto--;
	  
	  if (glen < 0)
	    glen = ((gto-gfrom+1 <= 6000) ? gto-gfrom+1 : 6000);
	  
	  InitPNG(PAR.getI("Output.resx"),   PAR.getI("Output.resy"),
		  PAR.getI("Output.offset"), gfrom, gto,
		  PAR.getI("Output.golap"), glen, grname);
	}
	
	MS->InitSensors(TheSeq);
    
	pred = Predict(TheSeq, MS);
    
	fprintf(stderr,"Optimal path length = %.4f\n",- pred->OptimalPath);
	
	if (graph)
	  pred->plotPred();
    
	Output(TheSeq, MS, pred, sequence, argc, argv, stdout);
	MS->PostAnalyse(pred);
    
	if (graph) {
	  fprintf(stderr,"\nDumping images (\"%s.---.png\")...", grname);
	  fflush(stderr);
	  ClosePNG();
	  fprintf(stderr, "done\n");
	}
	
	// Free used memory
	delete TheSeq; 
	delete MS;
	delete pred;//->resetPred();
	
      } // fin de traitement de chaque séquence....
  
      fprintf(stderr,"-------------------------------------");
      fprintf(stderr,"--------------------------------\n");
  
      return  0;
    }
}
