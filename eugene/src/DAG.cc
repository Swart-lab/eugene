// ------------------------------------------------------------------
// Copyright (C) 2004 INRA <eugene@ossau.toulouse.inra.fr>
//
// This program is open source; you can redistribute it and/or modify
// it under the terms of the Artistic License (see LICENSE file).
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
//
// You should have received a copy of Artistic License along with
// this program; if not, please see http://www.opensource.org
//
// $Id$
// ------------------------------------------------------------------
// File:     DAG.cc
// Contents: Modeling genes as a DAG and exploring it
// ------------------------------------------------------------------

#include "DAG.h"
#include <stdio.h>

extern Parameters PAR;

int verbose = 0;
int verbose2 = 0;

// ----------------------------------------------------------------
// Default constructor.
// ----------------------------------------------------------------
DAG :: DAG ()
{
  int i;

  StartPosition=INITIALSHIFT;
  EndPosition=INITIALSHIFT;
  EvidenceName[0]=0;
  pred = NULL;
  for (i = 0; i < NbTracks; i++) LBP[i].Path.SetState(i);
  NormalizingPath= 0.0;
  TheSeq = NULL;
}

// ----------------------------------------------------------------
// Default constructor.
// ----------------------------------------------------------------
DAG :: DAG (int start, int end, Parameters &PAR, DNASeq* Seq)
{
  StartPosition=start;
  EndPosition=end;
  strcpy(EvidenceName,"OPTIMAL");
  pred = NULL;
  NormalizingPath = 0.0;
  TheSeq = Seq;

  // Objectif -> limiter le nombre d'appel à la map de PAR
  SplicedStopPen = -PAR.getD("EuGene.SplicedStopPen");
  ExPrior    = PAR.getD("EuGene.ExonPrior");
  InPrior    = PAR.getD("EuGene.IntronPrior");
  IGPrior    = PAR.getD("EuGene.InterPrior"); 
  FivePrior  = PAR.getD("EuGene.FivePrimePrior");
  ThreePrior = PAR.getD("EuGene.ThreePrimePrior");
  IntronFivePrior = InPrior;
  MinCDSLen = PAR.getI("Output.MinCDSLen");
  estuse = PAR.getI("Sensor.Est.use");
}

// ----------------------------------------------------------------
// Constructor from another DAG. 
// Plug the extremities of the new one on pos start and end on the other
// ----------------------------------------------------------------
DAG :: DAG (int start, int end, DAG *RefDag,char* name)
{
  TheSeq = RefDag->TheSeq;
   StartPosition=start;
   EndPosition=end;
   pred = NULL;
   strcpy(EvidenceName,name);
   NormalizingPath = RefDag->NormalizingPath;

   //   for (int i = 0; i < NbTracks; i++) {
   //    LBP[i] = new Track(RefDag->LBP[i]);
   //  }

   // Objectif -> limiter le nombre d'appel à la map de PAR
   SplicedStopPen = RefDag->SplicedStopPen;
   ExPrior    = RefDag->ExPrior   ;
   InPrior    = RefDag->InPrior   ;
   IGPrior    = RefDag->IGPrior   ;
   FivePrior  = RefDag->FivePrior ;
   ThreePrior = RefDag->ThreePrior;

   MinCDSLen = RefDag->MinCDSLen;
   estuse =   RefDag->estuse ;
}

// ----------------------------------------------------------------
//  Default destructor.
// ----------------------------------------------------------------
DAG :: ~DAG ()
{
  for  (int i = 0;  i < NbTracks;  i ++) LBP[i].Zap();
}

// ---------------------------------------------------------------------------
// Initial and terminal costs
// ---------------------------------------------------------------------------
void DAG :: WeightThePrior()
{
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
  
  LBP[IntronF2T].Update(log(InPrior/6.0)/2.0);
  LBP[IntronF3TG].Update(log(InPrior/6.0)/2.0);
  LBP[IntronF3TA].Update(log(InPrior/6.0)/2.0);
  LBP[IntronR2AG].Update(log(InPrior/6.0)/2.0);
  LBP[IntronR3G].Update(log(InPrior/6.0)/2.0);
  LBP[IntronR3A].Update(log(InPrior/6.0)/2.0);

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
}

// ---------------------------------------------------------------------------
// Cost distributions on lengths
// ---------------------------------------------------------------------------
void DAG :: LoadDistLength()
{
  int j;

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
  
  LBP[IntronF2T].LoadPenalty(PAR.getC("EuGene.IntronDist"));
  LBP[IntronF3TG].LoadPenalty(PAR.getC("EuGene.IntronDist"));
  LBP[IntronF3TA].LoadPenalty(PAR.getC("EuGene.IntronDist"));
  LBP[IntronR2AG].LoadPenalty(PAR.getC("EuGene.IntronDist"));
  LBP[IntronR3G].LoadPenalty(PAR.getC("EuGene.IntronDist"));
  LBP[IntronR3A].LoadPenalty(PAR.getC("EuGene.IntronDist"));

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
}

// ----------------------------------------------------------------
//  Build the prediction by backtracing.
// ----------------------------------------------------------------
double DAG :: BuildPrediction (int Forward)
{
  double maxi,BestU,PBest[NbTracks];
  BackPoint *PrevBP[NbTracks];
  
  int j,k;

  // Insert best possible backpoint at the end of the algo (the
  // last insert is not automatically possible, cf cost dist. 
  // on length)
  for (j = 0; j < NbTracks; j++) {
    PrevBP[j] = LBP[j].BestUsable(TheSeq->SeqLen+1,&PBest[j],Forward,0);
    LBP[j].ForceNew(j,TheSeq->SeqLen+1,PBest[j],PrevBP[j]);
  }
  
  // Select where Backtrace should start
  j = 0;
  maxi = PBest[0];
  
  for (k = 1; k < NbTracks ; k++) {
    BestU = PBest[k];
    // Un test tordu pour casser le cou aux NaN
    if (isnan(maxi) || (BestU > maxi)) {
      maxi = BestU;
      j = k;
    }
  }

  pred = LBP[j].BackTrace(MinCDSLen,Forward);
  pred->OptimalPath = maxi+NormalizingPath;
  return maxi+NormalizingPath;
}

// ----------------------------------------------------------------
//  Build a prediction by reverting another one
// ----------------------------------------------------------------
void DAG :: BuildReversePrediction (DAG* dag)
{
  this->pred=new Prediction();
  for (int i= (dag->pred->size())-1; i>=0; i--)
    this->pred->add(dag->pred->getPos(i),dag->pred->getState(i));
}

// ----------------------------------------------------------------
//  Print infos of the DAG
// ----------------------------------------------------------------
void DAG :: Print()
{
  //  fprintf(stdout,"DAG %s normpath=%f start=%d end=%d HKey=%f CKey=%d\n",
  //	  EvidenceName,NormalizingPath,StartPosition,EndPosition,pred->getHashKey(),pred->getcodingdiffkey());
  //  fflush(stdout);
}

// ----------------------------------------------------------------
//  Shortest Path Algorithm with length constraint
// ----------------------------------------------------------------
void DAG :: ShortestPathAlgoForward (int position, DATA Data)
{
  
  int k;
  int Data_Len = TheSeq->SeqLen;
  double BestU, maxi = -NINFINITY;
  double PBest[NbTracks];
  BackPoint *PrevBP[NbTracks];
  signed   char best = 'Z';

  // Get information on possible spliced stops
  int StopStop = TheSeq->IsStopStop(position);
  int StartStop = TheSeq->IsStartStop(position);
    
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
#define INEED(K) if (!PrevBP[K]) PrevBP[K] = LBP[K].BestUsable(position, &PBest[K])
#define PICOMP(C,S,B,O) if (C && ISPOSSIBLE(S,B)) {\
	BestU = PBest[O]+Data.sig[DATA::S].weight[Signal::B];\
	if (isnan(maxi) || (BestU > maxi)) {maxi = BestU; best = O;}}
#define PICOMPEN(C,S,B,O,P) if (C && ISPOSSIBLE(S,B)) {\
      BestU = PBest[O]+Data.sig[DATA::S].weight[Signal::B]+(P);		\
	if (isnan(maxi) || (BestU > maxi)) {maxi = BestU; best = O;}}
#define INSERTANDPAYTHESLOPE(P,C)					\
    LBP[P].InsertNew(best, position, maxi,PrevBP[best]);\
    if (position < Data_Len){\
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

    // Exons F (No splice)
    for (k = 0; k <= 3; k++) {
      if ((position % 3 == k) && ISPOSSIBLE(Start,Forward))             INEED(UTR5F);
      if (ISPOSSIBLE(Ins,Forward)) { 	                           INEED(InitF1+(k+2)%3);
								   INEED(IntrF1+(k+2)%3);
      								   INEED(TermF1+(k+2)%3); }
      if (ISPOSSIBLE(Del,Forward)) {	                           INEED(InitF1+(k+1)%3);
      								   INEED(IntrF1+(k+1)%3);
      								   INEED(TermF1+(k+1)%3); }
    }

    // Exons F  (splicing, with and without spliced stops)
    if (ISPOSSIBLE(Acc,Forward)) {
      //      if (!(StopStop & DNASeq::StopAfterT))
	INEED(IntronF2T);
      
	//      if (!(StopStop & DNASeq::StopAfterTG))
	INEED(IntronF3TG);
   
	//      if (!(StopStop & DNASeq::StopAfterTA)) 
	INEED(IntronF3TA);
      
      // No spliced stops
      INEED(IntronF1+(position%3));
      INEED(IntronF1+((position+1)%3));
      INEED(IntronF1+((position+2)%3));
    }

    // Exons R  (No Splice)
    for (k = 3; k <= 6; k++) {
      if (((Data_Len-position) % 3 == k-3) && ISPOSSIBLE(Stop,Reverse)) INEED(UTR3R);
      if (ISPOSSIBLE(Ins,Reverse)) {				   INEED(InitR1+(k+2)%3);
								   INEED(IntrR1+(k+2)%3);
							  	   INEED(TermR1+(k+2)%3); }
      if (ISPOSSIBLE(Del,Reverse)) {				   INEED(InitR1+(k+1)%3);
								   INEED(IntrR1+(k+1)%3);
								   INEED(TermR1+(k+1)%3); }
    }
  
    // Exons R (splicing, with and without spliced stops)
    if (ISPOSSIBLE(Don,Reverse)) {
      //      if (!(StartStop & DNASeq::StopAfterG))
	INEED(IntronR3G);
      
	//      if (!(StartStop & DNASeq::StopAfterA))
	INEED(IntronR3A);
   
	//      if (!(StartStop & DNASeq::StopAfterAG)) 
	INEED(IntronR2AG);

      // No spliced stops
      INEED(IntronR1+((Data_Len-position) % 3));
      INEED(IntronR1+((Data_Len-position+2) % 3));
      INEED(IntronR1+((Data_Len-position+1) % 3));
    }

    // Introns F et R (with no spliceable stops)
    if (ISPOSSIBLE(Don,Forward)) 
      for (k=0; k<3; k++) {	
	INEED(InitF1+((position-k+3) % 3));
	INEED(IntrF1+((position-k+3) % 3)); 
      }

    if (ISPOSSIBLE(Acc,Reverse)) 
      for (k=0; k<3; k++) {
      INEED(IntrR1+((Data_Len-position-k+3) % 3));
      INEED(TermR1+((Data_Len-position-k+3) % 3));
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
    if (ISPOSSIBLE(Stop,Forward)) {			           INEED(TermF1+position%3);
			   			                   INEED(SnglF1+position%3); }
    // UTR 5'reverse
    if (ISPOSSIBLE(Start,Reverse)) {                               INEED(InitR1+((Data_Len-position) % 3));
    				                                   INEED(SnglR1+((Data_Len-position) % 3)); }
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
      PICOMP((position % 3 == k),Start,Forward,UTR5F);
      // Il y a une insertion (frameshift). Saut de positionléotide ignore.
      PICOMP(true,Ins,Forward,InitF1+(k+2)%3);
      // Il y a une deletion (frameshift)
      PICOMP(true,Del,Forward,InitF1+(k+1)%3);
    
      // On va tout droit.
      // S'il y  a un STOP en phase on ne peut continuer
      if ((position % 3 == k))
	LBP[InitF1+k].Update(Data.sig[DATA::Stop].weight[Signal::ForwardNo]);
      LBP[InitF1+k].Update(Data.sig[DATA::Don].weight[Signal::ForwardNo]);

      INSERTANDPAYTHESLOPE(InitF1+k,k);
    }
    // ----------------------------------------------------------------
    // ------------------------- Inits en reverse ---------------------
    // ----------------------------------------------------------------
    for (k = 3; k<6; k++) {
      maxi = NINFINITY;
      
      // - on recommence a coder (Donneur) Ca vient d'un intron (no spliceable stop)
      PICOMPEN(true,Don,Reverse,IntronR1+((Data_Len-position-k) % 3),((((Data_Len-position) % 3) == k-3) ? Data.sig[DATA::Stop].weight[Signal::ReverseNo] : 0.0));
      // Not AfterG
      PICOMPEN(((Data_Len-position-k) % 3) == 2,Don,Reverse,IntronR3G ,((StartStop & DNASeq::StopAfterG)  ? SplicedStopPen : 0.0));
      // Not AfterA
      PICOMPEN(((Data_Len-position-k) % 3) == 2,Don,Reverse,IntronR3A ,((StartStop & DNASeq::StopAfterA)  ? SplicedStopPen : 0.0));
      // Not AfterAG
      PICOMPEN(((Data_Len-position-k) % 3) == 1,Don,Reverse,IntronR2AG,((StartStop & DNASeq::StopAfterAG) ? SplicedStopPen : 0.0));

      // Il y a une insertion (frameshift)
      PICOMP(true,Ins,Reverse, InitR1+(k+2)%3);
      // Il y a une deletion (frameshift)
      PICOMP(true,Del,Reverse, InitR1+(k+1)%3);

      if (((Data_Len-position) % 3 == k-3)) 
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
      PICOMP((position % 3 == k),Start,Forward,UTR5F);
      // Il y a une insertion (frameshift). Saut de positionléotide ignore.
      PICOMP(true,Ins,Forward,SnglF1+(k+2)%3);
      // Il y a une deletion (frameshift)
      PICOMP(true,Del,Forward,SnglF1+(k+1)%3);
    
      // On va tout droit.
      // S'il y  a un STOP en phase on ne peut continuer
      if ((position % 3 == k))
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
      PICOMP(((Data_Len-position) % 3 == k-3),Stop,Reverse,UTR3R);
      // Il y a une insertion (frameshift)
      PICOMP(true,Ins,Reverse, SnglR1+(k+2)%3);
      // Il y a une deletion (frameshift)
      PICOMP(true,Del,Reverse, SnglR1+(k+1)%3);

      if (((Data_Len-position) % 3 == k-3)) 
	LBP[SnglF1+k].Update(Data.sig[DATA::Stop].weight[Signal::ReverseNo]);
      LBP[SnglF1+k].Update(Data.sig[DATA::Don].weight[Signal::ReverseNo]);

      INSERTANDPAYTHESLOPE(SnglF1+k,k);
    }
    // ----------------------------------------------------------------
    // ------------------ Intrs en forward ----------------------------
    // ----------------------------------------------------------------
    for (k = 0; k < 3; k++) {
      maxi = NINFINITY;     

      // On recommence a coder (Accepteur). Ca vient d'un intron (no spliceable stop)
      PICOMP(true,Acc,Forward,IntronF1+((position-k+3) % 3));
      // Not AfterT
      PICOMPEN(((position-k+3) % 3) == 1,Acc,Forward,IntronF2T ,((StopStop & DNASeq::StopAfterT)  ? SplicedStopPen : 0.0));
      // Not AfterTG
      PICOMPEN(((position-k+3) % 3) == 2,Acc,Forward,IntronF3TG,((StopStop & DNASeq::StopAfterTG) ? SplicedStopPen : 0.0));
      // Not AfterTA
      PICOMPEN(((position-k+3) % 3) == 2,Acc,Forward,IntronF3TA,((StopStop & DNASeq::StopAfterTA) ? SplicedStopPen : 0.0));
      // Il y a une insertion (frameshift). Saut de positionléotide ignore.
      PICOMP(true,Ins,Forward,IntrF1+(k+2)%3);
      // Il y a une deletion (frameshift)
      PICOMP(true,Del,Forward,IntrF1+(k+1)%3);
    
      // On va tout droit.
      // S'il y  a un STOP en phase on ne peut continuer
      if ((position % 3 == k))
	LBP[IntrF1+k].Update(Data.sig[DATA::Stop].weight[Signal::ForwardNo]);
      LBP[IntrF1+k].Update(Data.sig[DATA::Don].weight[Signal::ForwardNo]);

      INSERTANDPAYTHESLOPE(IntrF1+k,k);
    }
    // ----------------------------------------------------------------
    // ------------------------- Intrs en reverse ---------------------
    // ----------------------------------------------------------------
    for (k = 3; k<6; k++) {
      maxi = NINFINITY;
      
      // - on recommence a coder (Donneur) Ca vient d'un intron (no spliceable stop)
      PICOMPEN(true,Don,Reverse,IntronR1+((Data_Len-position-k) % 3),((((Data_Len-position) % 3) == k-3) ? Data.sig[DATA::Stop].weight[Signal::ReverseNo] : 0.0));
      // Not AfterG
      PICOMPEN(((Data_Len-position-k) % 3) == 2,Don,Reverse,IntronR3G ,((StartStop & DNASeq::StopAfterG)  ? SplicedStopPen : 0.0));
      // Not AfterA
      PICOMPEN(((Data_Len-position-k) % 3) == 2,Don,Reverse,IntronR3A ,((StartStop & DNASeq::StopAfterA)  ? SplicedStopPen : 0.0));
      // Not AfterAG
      PICOMPEN(((Data_Len-position-k) % 3) == 1,Don,Reverse,IntronR2AG,((StartStop & DNASeq::StopAfterAG) ? SplicedStopPen : 0.0));
        // Il y a une insertion (frameshift)
      PICOMP(true,Ins,Reverse, IntrR1+(k+2)%3);
      // Il y a une deletion (frameshift)
      PICOMP(true,Del,Reverse, IntrR1+(k+1)%3);

      if (((Data_Len-position) % 3 == k-3)) 
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
      PICOMP(true,Acc,Forward,IntronF1+((position-k+3) % 3));
      // Not AfterT
      PICOMPEN(((position-k+3) % 3) == 1,Acc,Forward,IntronF2T ,((StopStop & DNASeq::StopAfterT)  ? SplicedStopPen :0.0));
      // Not AfterTG
      PICOMPEN(((position-k+3) % 3) == 2,Acc,Forward,IntronF3TG,((StopStop & DNASeq::StopAfterTG) ? SplicedStopPen : 0.0));
      // Not AfterTA
      PICOMPEN(((position-k+3) % 3) == 2,Acc,Forward,IntronF3TA,((StopStop & DNASeq::StopAfterTA) ? SplicedStopPen : 0.0));
      // Il y a une insertion (frameshift). Saut de positionléotide ignore.
      PICOMP(true,Ins,Forward,TermF1+(k+2)%3);
      // Il y a une deletion (frameshift)
      PICOMP(true,Del,Forward,TermF1+(k+1)%3);
    
      // On va tout droit.
      // S'il y  a un STOP en phase on ne peut continuer
      if ((position % 3 == k))
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
      PICOMP(((Data_Len-position) % 3 == k-3),Stop,Reverse,UTR3R);
      // Il y a une insertion (frameshift)
      PICOMP(true,Ins,Reverse, TermR1+(k+2)%3);
      // Il y a une deletion (frameshift)
      PICOMP(true,Del,Reverse, TermR1+(k+1)%3);

      if (((Data_Len-position) % 3 == k-3)) 
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
    LBP[InterGen].InsertNew(best, position, maxi,PrevBP[best]);

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
    PICOMP(true,Stop,Forward, TermF1+position%3);
    PICOMP(true,Stop,Forward, SnglF1+position%3);
    // On peut venir aussi d'un intron d'UTR.
    PICOMP(true,Acc,Forward, IntronU3F);

    // On reste 3' direct
    LBP[UTR3F].InsertNew(best, position, maxi,PrevBP[best]);

    INSERTANDPAYTHESLOPE(UTR3F,DATA::UTR3F);
    // ----------------------------------------------------------------
    // ----------------------- UTR 5'reverse --------------------------
    // ----------------------------------------------------------------
    maxi = NINFINITY;
    
    // Ca vient d'un Sngl ou Init reverse + START
    PICOMP(true,Start,Reverse, InitR1+((Data_Len-position) % 3));
    PICOMP(true,Start,Reverse, SnglR1+((Data_Len-position) % 3));
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
    LBP[UTR3R].InsertNew(best, position, maxi,PrevBP[best]);

    INSERTANDPAYTHESLOPE(UTR3R,DATA::UTR3R);
    // ----------------------------------------------------------------
    // ---------------- Introns de phase k forward --------------------
    // ----------------------------------------------------------------
    for (k = 0; k<3; k++) {
      maxi = NINFINITY;
      
      // - on quitte un Init ou un Intr
      // no spliceable stop: 
      if (!(((StartStop & DNASeq::StopAfterT) && k == 1) ||
			((StartStop & (DNASeq::StopAfterTG|DNASeq::StopAfterTA)) && k == 2)))
	  {
		  PICOMPEN(true,Don,Forward, InitF1+((position-k+3) % 3),(k == 0 ? Data.sig[DATA::Stop].weight[Signal::ForwardNo] : 0.0));
		  PICOMPEN(true,Don,Forward, IntrF1+((position-k+3) % 3),(k == 0 ? Data.sig[DATA::Stop].weight[Signal::ForwardNo] : 0.0));
	  }
      // On reste intronique
      LBP[IntronF1+k].Update(Data.sig[DATA::Acc].weight[Signal::ForwardNo]);
      
      INSERTANDPAYTHESLOPE(IntronF1+k,DATA::IntronF);
    }
    // ----------------------------------------------------------------
    // ----------------- Introns forward speciaux ---------------------
    // ----------------------------------------------------------------
    //
    // --- Intron Phase 1 after a T (GA|AA|AG)
    //
    maxi = NINFINITY;
    k = 1;
    if (StartStop & DNASeq::StopAfterT) {
      // - on quitte un Init ou un Intr
      PICOMP(true,Don,Forward, InitF1+((position-k+3) % 3));
      PICOMP(true,Don,Forward, IntrF1+((position-k+3) % 3));
    }
    // On reste intronique
    LBP[IntronF2T].Update(Data.sig[DATA::Acc].weight[Signal::ForwardNo]);
    INSERTANDPAYTHESLOPE(IntronF2T,DATA::IntronF);

    //
    // --- Intron Phase 2 after an TG|A
    //
    maxi = NINFINITY;
    k = 2;
    if (StartStop & DNASeq::StopAfterTG) {
      // - on quitte un Init ou un Intr
      PICOMP(true,Don,Forward, InitF1+((position-k+3) % 3));
      PICOMP(true,Don,Forward, IntrF1+((position-k+3) % 3));
    }
    // On reste intronique
    LBP[IntronF3TG].Update(Data.sig[DATA::Acc].weight[Signal::ForwardNo]);
    INSERTANDPAYTHESLOPE(IntronF3TG,DATA::IntronF);

    //
    // --- Intron Phase 2 after a TA(A|G)
    //
    maxi = NINFINITY;
    k = 2;
    if (StartStop & DNASeq::StopAfterTA) {
      // - on quitte un Init ou un Intr
      PICOMP(true,Don,Forward, InitF1+((position-k+3) % 3));
      PICOMP(true,Don,Forward, IntrF1+((position-k+3) % 3));
    }
    // On reste intronique
    LBP[IntronF3TA].Update(Data.sig[DATA::Acc].weight[Signal::ForwardNo]);
    INSERTANDPAYTHESLOPE(IntronF3TA,DATA::IntronF);

    // ----------------------------------------------------------------
    // ----------------- Introns de phase -k reverse ------------------
    // ----------------------------------------------------------------
    for (k = 0; k<3; k++) {
      maxi = NINFINITY;
      
      // On quitte un Intr ou un Term reverse
      // no spliceable stop: 
      if (!(((StopStop & (DNASeq::StopAfterG|DNASeq::StopAfterA)) && k == 2) ||
	    ((StopStop & DNASeq::StopAfterAG) && k == 1)))  {
	PICOMP(true,Acc,Reverse, IntrR1+((Data_Len-position-k) % 3));
	PICOMP(true,Acc,Reverse, TermR1+((Data_Len-position-k) % 3));
      }
      // On reste intronique
      LBP[IntronR1+k].Update(Data.sig[DATA::Acc].weight[Signal::ReverseNo]);

      INSERTANDPAYTHESLOPE(IntronR1+k,DATA::IntronR);
    }
    // ----------------------------------------------------------------
    // ----------------- Introns reverse speciaux ---------------------
    // ----------------------------------------------------------------
    //
    // --- Intron Phase 1 after a G (AT)
    //
    maxi = NINFINITY;
    k = 2;
    if (StopStop & DNASeq::StopAfterG) {
      // - on quitte un Intr ou un Term
      PICOMP(true,Acc,Reverse, IntrR1+((Data_Len-position-k) % 3));
      PICOMP(true,Acc,Reverse, TermR1+((Data_Len-position-k) % 3));
    }
    // On reste intronique
    LBP[IntronR3G].Update(Data.sig[DATA::Acc].weight[Signal::ReverseNo]);
    INSERTANDPAYTHESLOPE(IntronR3G,DATA::IntronR);

    //
    // --- Intron Phase 1 after an A (GT|AT)
    //
    maxi = NINFINITY;
    k = 2;
    if (StopStop & DNASeq::StopAfterA) {
      // - on quitte un Intr ou un Term
      PICOMP(true,Acc,Reverse, IntrR1+((Data_Len-position-k) % 3));
      PICOMP(true,Acc,Reverse, TermR1+((Data_Len-position-k) % 3));
    }
    // On reste intronique
    LBP[IntronR3A].Update(Data.sig[DATA::Acc].weight[Signal::ReverseNo]);
    INSERTANDPAYTHESLOPE(IntronR3A,DATA::IntronR);

    //
    // --- Intron Phase 2 after an AG, AA ou GA (T)
    //
    maxi = NINFINITY;
    k = 1;
    if (StopStop & DNASeq::StopAfterAG) {
      // - on quitte un Intr ou un Term
      PICOMP(true,Acc,Reverse, IntrR1+((Data_Len-position-k) % 3));
      PICOMP(true,Acc,Reverse, TermR1+((Data_Len-position-k) % 3));
    }
    // On reste intronique
    LBP[IntronR2AG].Update(Data.sig[DATA::Acc].weight[Signal::ReverseNo]);
    INSERTANDPAYTHESLOPE(IntronR2AG,DATA::IntronR);
}


// ----------------------------------------------------------------
//  Shortest Path Algorithm with length constraint (reverse)
// ----------------------------------------------------------------
void DAG :: ShortestPathAlgoBackward (int position, DATA Data, int NoContentsUpdate)
{

  int k;
  int Data_Len = TheSeq->SeqLen;
  double BestU, maxi = -NINFINITY;
  double PBest[NbTracks];
  BackPoint *PrevBP[NbTracks];
  signed   char best = 'Z';

  // Get information on possible spliced stops
  int StopStop = TheSeq->IsStopStop(position);
  int StartStop = TheSeq->IsStartStop(position);
    
  // normalisation des couts 
  maxi = -NINFINITY;
  for (k = 0 ; k < NbTracks; k++) {
    BestU = LBP[k].Optimal;
    if ((BestU > NINFINITY) && (BestU < maxi)) maxi =BestU;
  }

  for (k = 0 ; k < NbTracks; k++) 
    LBP[k].Update(-maxi);
  NormalizingPath += maxi;

#define INEEDR(K) if (!PrevBP[K]) PrevBP[K] = LBP[K].BestUsable(position, &PBest[K],1)
#define INSERTANDPAYTHESLOPER(P,C)					\
    LBP[P].InsertNew(best, position, maxi,PrevBP[best]);\
    if (position < Data_Len){\
      if (!NoContentsUpdate) {\
	LBP[P].Update(Data.contents[C]);	\
	LBP[P].PayTheSlope();			\
      }\
    }

    // ----------------------------------------------------------------
    // Calcul des meilleures opening edges
    // ----------------------------------------------------------------
    for (k = 0; k < NbTracks; k++) {
      PrevBP[k] = NULL;
      PBest[k] = NINFINITY;
    }
    // ---------------------------
    // 1- Exons F spliced (from Intron)
    // ---------------------------
    if (ISPOSSIBLE(Don,Forward)) {
      INEEDR(IntronF1);
      INEEDR(IntronF2);
      INEEDR(IntronF3);
      INEEDR(IntronF2T);
      INEEDR(IntronF3TG);
      INEEDR(IntronF3TA);
    }
    // ---------------------------
    // 2- Exons F (from UTR3F)
    // ---------------------------
    if (ISPOSSIBLE(Stop,Forward))           
      INEEDR(UTR3F);
    // ---------------------------
    // 3- Intron F (from ExonF)
    // ---------------------------
    if (ISPOSSIBLE(Acc,Forward)) {
      INEEDR(IntrF1);
      INEEDR(IntrF2);
      INEEDR(IntrF3);
      INEEDR(TermF1);
      INEEDR(TermF2);
      INEEDR(TermF3);
    }
    // ---------------------------
    // 4- UTR5F (from SnglF) - check
    // ---------------------------
    if (ISPOSSIBLE(Start,Forward)) {
      INEEDR(SnglF1+((position+1)%3));
      INEEDR(InitF1+((position+1)%3));
    }
    // ---------------------------
    // 5- UTR5F (from IntronU5F)
    // ---------------------------
    if (ISPOSSIBLE(Don,Forward))
      INEEDR(IntronU5F);
    // ---------------------------
    // 6- IntronU5F (from UTR5F)
    // ---------------------------
    if (ISPOSSIBLE(Acc,Forward))
      INEEDR(UTR5F);
    // ---------------------------
    // 7- UTR3F (from InterGen)
    // ---------------------------
    if (ISPOSSIBLE(tStart,Forward))
      INEEDR(InterGen);
    // ---------------------------
    // 8- UTR3F (from IntronU3F)
    // ---------------------------
    if (ISPOSSIBLE(Don,Forward))
      INEEDR(IntronU3F);
    // ---------------------------
    // 9- IntronU3F (from UTR3F)
    // ---------------------------
    if (ISPOSSIBLE(Acc,Forward))
      INEEDR(UTR3F);
    // ---------------------------
    // 10- InterG (from UTR5F)
    // ---------------------------
    if (ISPOSSIBLE(tStart,Forward))
      INEEDR(UTR5F);
    // ---------------------------
    // 11- InterG (from UTR3R)
    // ---------------------------
    if (ISPOSSIBLE(tStop,Reverse))
      INEEDR(UTR3R);
    // ---------------------------
    // 12- IntronU5R (from UTR5R)
    // ---------------------------
    if (ISPOSSIBLE(Don,Reverse))
      INEEDR(UTR5R);
    // ---------------------------
    // 13- UTR5R (from InterGen)
    // ---------------------------
    if (ISPOSSIBLE(tStart,Reverse))
      INEEDR(InterGen);
    // ---------------------------
    // 14- UTR5R (from IntronU5R)
    // ---------------------------
    if (ISPOSSIBLE(Acc,Reverse))
      INEEDR(IntronU5R);
    // ---------------------------
    // 15- IntronU3R (from UTR3R)
    // ---------------------------
    if (ISPOSSIBLE(Don,Reverse))
      INEEDR(UTR3R);
    // ---------------------------
    // 16- UTR3R (from IntronU3R)
    // ---------------------------
    if (ISPOSSIBLE(Acc,Reverse))
      INEEDR(IntronU3R);
    // ---------------------------
    // 17- UTR3R (from ExonR)
    // ---------------------------
    if (ISPOSSIBLE(Stop,Reverse)) {
      INEEDR(SnglR1+(Data_Len-position-1)%3);
      INEEDR(TermR1+(Data_Len-position-1)%3);
    }
    // ---------------------------
    // 18- ExonR (from UTR5R)
    // ---------------------------
    if (ISPOSSIBLE(Start,Reverse))
      INEEDR(UTR5R);
    // ---------------------------
    // 19- ExonR (from IntronR)
    // ---------------------------
    if (ISPOSSIBLE(Acc,Reverse)){
      INEEDR(IntronR1);
      INEEDR(IntronR2);
      INEEDR(IntronR3);
      INEEDR(IntronR2AG);
      INEEDR(IntronR3A);
      INEEDR(IntronR3G);
    }
    // ---------------------------
    // 20- IntronR (from ExonR)
    // ---------------------------
    if (ISPOSSIBLE(Don,Reverse)) {
      INEEDR(InitR1);
      INEEDR(InitR2);
      INEEDR(InitR3);
      INEEDR(IntrR1);
      INEEDR(IntrR2);
      INEEDR(IntrR3);
    }

    // ----------------------------------------------------------------
    // ------------------ Inits en forward ----------------------------
    // ----------------------------------------------------------------
    for (k = 0; k < 3; k++) {
      maxi = NINFINITY;     
      
      // Début à partit d'un IntronF
      PICOMP(true,Don,Forward,IntronF1+(position-k+4)%3);
      // TODO 1- spliceable STOPs
      // TODO 2- prise en compte STOP dans le PICOMP (cf. Init Reverse en AlgoForward)
    
      // On va tout droit.
      // S'il y  a un STOP en phase on ne peut continuer
      if (((position+1)%3 == k))
	LBP[InitF1+k].Update(Data.sig[DATA::Stop].weight[Signal::ForwardNo]);
      LBP[InitF1+k].Update(Data.sig[DATA::Don].weight[Signal::ForwardNo]);

      INSERTANDPAYTHESLOPER(InitF1+k,k);
    }
    // ----------------------------------------------------------------
    // ------------------------- Inits en reverse ---------------------
    // ----------------------------------------------------------------
    for (k = 0; k<3; k++) {
      maxi = NINFINITY;
      
      // - on recommence a coder (Start) d'une UTR5R
      PICOMP((Data_Len-position-1)%3 == k,Start,Reverse,UTR5R);
      // TODO prise en compte STOP dans le PICOMP (cf Init Reverse en Forward)

      if (((Data_Len-position-1) % 3 == k)) 
	LBP[InitR1+k].Update(Data.sig[DATA::Stop].weight[Signal::ReverseNo]);
      LBP[InitR1+k].Update(Data.sig[DATA::Don].weight[Signal::ReverseNo]);

      INSERTANDPAYTHESLOPE(InitR1+k,k+3);
    }
    // ----------------------------------------------------------------
    // ------------------ Sngl en forward ----------------------------
    // ----------------------------------------------------------------
    for (k = 0; k < 3; k++) {
      maxi = NINFINITY;     
      
      // On commence a coder (Stop) sur une UTR3F
      PICOMP(((position+1)%3 == k),Stop,Forward,UTR3F);
      // TODO prise en compte STOP dans le PICOMP (cf Init Reverse en Forward)
      
      // On va tout droit.
      // S'il y  a un STOP en phase on ne peut continuer
      if (((position+1)%3 == k))
	LBP[SnglF1+k].Update(Data.sig[DATA::Stop].weight[Signal::ForwardNo]);
      LBP[SnglF1+k].Update(Data.sig[DATA::Don].weight[Signal::ForwardNo]);
      
      INSERTANDPAYTHESLOPE(SnglF1+k,k);
    }
    // ----------------------------------------------------------------
    // ------------------------- Sngl en reverse ---------------------
    // ----------------------------------------------------------------
    for (k = 0; k<3; k++) {
      maxi = NINFINITY;
      
      // On commence a coder (Start). Ca vient d'une UTR 5' reverse
      PICOMP(((Data_Len-position-1)%3 == k),Start,Reverse,UTR5R);
      // TODO prise en compte STOP dans le PICOMP (cf Init Reverse en Forward)

      if (((Data_Len-position-1)%3 == k)) 
	LBP[SnglR1+k].Update(Data.sig[DATA::Stop].weight[Signal::ReverseNo]);
      LBP[SnglR1+k].Update(Data.sig[DATA::Don].weight[Signal::ReverseNo]);

      INSERTANDPAYTHESLOPE(SnglR1+k,k+3);
    }
    // ----------------------------------------------------------------
    // ------------------ Intrs en forward ----------------------------
    // ----------------------------------------------------------------
    for (k = 0; k < 3; k++) {
      maxi = NINFINITY;     

      // On recommence a coder (Donneur). Ca vient d'un intronF
      PICOMP(true,Don,Forward,IntronF1+((position-k+4) % 3));
      // TODO 1- spliceable STOPs
      // TODO 2- prise en compte STOP dans le PICOMP (cf. Init Reverse en AlgoForward)

      // On va tout droit.
      // S'il y  a un STOP en phase on ne peut continuer
      if (((position+1)%3 == k))
	LBP[IntrF1+k].Update(Data.sig[DATA::Stop].weight[Signal::ForwardNo]);
      LBP[IntrF1+k].Update(Data.sig[DATA::Don].weight[Signal::ForwardNo]);

      INSERTANDPAYTHESLOPE(IntrF1+k,k);
    }
    // ----------------------------------------------------------------
    // ------------------------- Intrs en reverse ---------------------
    // ----------------------------------------------------------------
    for (k = 0; k<3; k++) {
      maxi = NINFINITY;
      
      // - on recommence a coder (Acc) Ca vient d'un intronR
      PICOMP(true,Acc,Reverse,IntronR1+((Data_Len-position-k-1)%3));
      // TODO 1- spliceable STOPs
      // TODO 2- prise en compte STOP dans le PICOMP (cf. Init Reverse en AlgoForward)

      if (((Data_Len-position-1)%3 == k)) 
	LBP[IntrF1+k].Update(Data.sig[DATA::Stop].weight[Signal::ReverseNo]);
      LBP[IntrF1+k].Update(Data.sig[DATA::Don].weight[Signal::ReverseNo]);

      INSERTANDPAYTHESLOPE(IntrR1+k,k+3);
    }
    // ----------------------------------------------------------------
    // ------------------ Terms en forward ----------------------------
    // ----------------------------------------------------------------
    for (k = 0; k < 3; k++) {
      maxi = NINFINITY;     
      
      // On recommence a coder (Stop). Ca vient d'une UTR3R
      PICOMP(((Data_Len-position-1)%3 == k),Stop,Forward,UTR3R);
      // TODO 2- prise en compte STOP dans le PICOMP (cf. Init Reverse en AlgoForward)

      // On va tout droit.
      // S'il y  a un STOP en phase on ne peut continuer
      if (((position+1)%3 == k))
	LBP[TermF1+k].Update(Data.sig[DATA::Stop].weight[Signal::ForwardNo]);
      LBP[TermF1+k].Update(Data.sig[DATA::Don].weight[Signal::ForwardNo]);

      INSERTANDPAYTHESLOPE(TermF1+k,k);
    }
    // ----------------------------------------------------------------
    // ------------------------- Terms en reverse ---------------------
    // ----------------------------------------------------------------
    for (k = 0; k<3; k++) {
      maxi = NINFINITY;
      
      // On commence a coder (Acc). Ca vient d'un intronR
      PICOMP(true,Acc,Reverse,IntronR1+(Data_Len-position-1-k)%3);

      if (((Data_Len-position-1)%3 == k)) 
	LBP[TermR1+k].Update(Data.sig[DATA::Stop].weight[Signal::ReverseNo]);
      LBP[TermR1+k].Update(Data.sig[DATA::Don].weight[Signal::ReverseNo]);

      INSERTANDPAYTHESLOPE(TermR1+k,k+3);
    }
    // ----------------------------------------------------------------
    // ------------------------ Intergenique --------------------------
    // ----------------------------------------------------------------
    // Ca peut venir d'une fin de 3' direct ou de 5' reverse
    maxi = NINFINITY;
    
    // From 5' forward
    PICOMP(true,tStart,Forward, UTR5F);
    // From 3' reverse
    PICOMP(true,tStop,Reverse, UTR3R);
 
    // On reste intergenique
    // et les transstartNO/TransstopNO ???
    LBP[InterGen].InsertNew(best, position, maxi,PrevBP[best]);

    INSERTANDPAYTHESLOPE(InterGen,DATA::InterG);
    // ----------------------------------------------------------------
    // ---------------------- UTR 5' direct ---------------------------
    // ----------------------------------------------------------------
    maxi = NINFINITY;
    
    // On vient d'un exon (SnglF ou InitF) sur un Start
    PICOMP(true,Start,Forward, SnglF1+(position+1)%3);
    // On peut venir aussi d'un intron d'UTR.
    PICOMP(true,Don,Forward, IntronU5F);

    // On reste 5' direct. On ne prend pas le Start eventuel.
    LBP[UTR5F].Update(Data.sig[DATA::Start].weight[Signal::ForwardNo]);

    INSERTANDPAYTHESLOPE(UTR5F,DATA::UTR5F);
    // ----------------------------------------------------------------
    // ------------------- Introns d'UTR5F ----------------------------
    // ----------------------------------------------------------------
    maxi = NINFINITY;
      
    // on quitte l'UTR 
    PICOMP(true,Acc,Forward, UTR5F);

    // On reste intronique
    LBP[IntronU5F].Update(Data.sig[DATA::Acc].weight[Signal::ForwardNo]);

    INSERTANDPAYTHESLOPE(IntronU5F,DATA::IntronUTRF);
    // ----------------------------------------------------------------
    // ------------------- Introns d'UTR3F ----------------------------
    // ----------------------------------------------------------------
    maxi = NINFINITY;
      
    // on quitte l'UTR 
    PICOMP(true,Acc,Forward, UTR3F);

    // On reste intronique
    LBP[IntronU3F].Update(Data.sig[DATA::Acc].weight[Signal::ForwardNo]);

    INSERTANDPAYTHESLOPE(IntronU3F,DATA::IntronUTRF);
    // ----------------------------------------------------------------
    // ---------------------- UTR 3' direct ---------------------------
    // ----------------------------------------------------------------
    maxi = NINFINITY;

    // Ca vient d'InterG (tStop)
    PICOMP(true,tStop,Forward, InterGen);
    // On peut venir aussi d'un intron d'UTR.
    PICOMP(true,Don,Forward, IntronU3F);

    // On reste 3' direct
    LBP[UTR3F].InsertNew(best, position, maxi,PrevBP[best]);

    INSERTANDPAYTHESLOPE(UTR3F,DATA::UTR3F);
    // ----------------------------------------------------------------
    // ----------------------- UTR 5'reverse --------------------------
    // ----------------------------------------------------------------
    maxi = NINFINITY;
    
    // Ca vient d'interG (tStartR)
    PICOMP(true,tStart,Reverse, InterGen);
    // On peut venir aussi d'un intron d'UTR.
    PICOMP(true,Acc,Reverse, IntronU5R);

    // On reste 5' reverse
    LBP[UTR5R].Update(Data.sig[DATA::Start].weight[Signal::ReverseNo]);

    INSERTANDPAYTHESLOPE(UTR5R,DATA::UTR5R);
    // ----------------------------------------------------------------
    // ------------------- Introns d'UTR5R ----------------------------
    // ----------------------------------------------------------------
    maxi = NINFINITY;
    
    // On quitte une UTR5R
    PICOMP(true,Don,Reverse, UTR5R);

    // On reste intronique
    LBP[IntronU5R].Update(Data.sig[DATA::Acc].weight[Signal::ReverseNo]);

    INSERTANDPAYTHESLOPE(IntronU5R,DATA::IntronUTRR);
    // ----------------------------------------------------------------
    // ------------------- Introns d'UTR3R ----------------------------
    // ----------------------------------------------------------------
    maxi = NINFINITY;
    
    // On quitte une UTR5R
    PICOMP(true,Don,Reverse, UTR3R);

    // On reste intronique
    LBP[IntronU3R].Update(Data.sig[DATA::Acc].weight[Signal::ReverseNo]);

    INSERTANDPAYTHESLOPE(IntronU3R,DATA::IntronUTRR);
    // ----------------------------------------------------------------
    // ----------------------- UTR 3'reverse --------------------------
    // ----------------------------------------------------------------
    maxi = NINFINITY;
    
    // On demarre depuis Sngl ou Term (Stop)
    PICOMP(true,Stop,Reverse, SnglR1+(Data_Len-position-1)%3);
    PICOMP(true,Stop,Reverse, TermR1+(Data_Len-position-1)%3);
    // On peut venir aussi d'un intron d'UTR.
    PICOMP(true,Acc,Reverse, IntronU3R);

    // On reste 3' reverse
    LBP[UTR3R].InsertNew(best, position, maxi,PrevBP[best]);

    INSERTANDPAYTHESLOPE(UTR3R,DATA::UTR3R);
    // ----------------------------------------------------------------
    // ---------------- Introns de phase k forward --------------------
    // ----------------------------------------------------------------
    for (k = 0; k<3; k++) {
      maxi = NINFINITY;
      
      // - on quitte un Term ou un Intr sur un Acc
      // TODO: no spliceable stop
      PICOMP(true,Acc,Forward,IntrF1+(position-k+4)%3);
      PICOMP(true,Acc,Forward,TermF1+(position-k+4)%3);
    
      // On reste intronique
      LBP[IntronF1+k].Update(Data.sig[DATA::Acc].weight[Signal::ForwardNo]);
      
      INSERTANDPAYTHESLOPE(IntronF1+k,DATA::IntronF);
    }
    // ----------------------------------------------------------------
    // ----------------- Introns forward speciaux ---------------------
    // ----------------------------------------------------------------
    //
    // --- Intron Phase 1 after a T (GA|AA|AG)
    //
    /* TODO Spliceable stops
    maxi = NINFINITY;
    k = 1;
    if (StartStop & DNASeq::StopAfterT) {
      // - on quitte un Init ou un Intr
      PICOMP(true,Don,Forward, InitF1+((position-k+3) % 3));
      PICOMP(true,Don,Forward, IntrF1+((position-k+3) % 3));
    }
    // On reste intronique
    LBP[IntronF2T].Update(Data.sig[DATA::Acc].weight[Signal::ForwardNo]);
    INSERTANDPAYTHESLOPE(IntronF2T,DATA::IntronF);

    //
    // --- Intron Phase 2 after an TG|A
    //
    maxi = NINFINITY;
    k = 2;
    if (StartStop & DNASeq::StopAfterTG) {
      // - on quitte un Init ou un Intr
      PICOMP(true,Don,Forward, InitF1+((position-k+3) % 3));
      PICOMP(true,Don,Forward, IntrF1+((position-k+3) % 3));
    }
    // On reste intronique
    LBP[IntronF3TG].Update(Data.sig[DATA::Acc].weight[Signal::ForwardNo]);
    INSERTANDPAYTHESLOPE(IntronF3TG,DATA::IntronF);

    //
    // --- Intron Phase 2 after a TA(A|G)
    //
    maxi = NINFINITY;
    k = 2;
    if (StartStop & DNASeq::StopAfterTA) {
      // - on quitte un Init ou un Intr
      PICOMP(true,Don,Forward, InitF1+((position-k+3) % 3));
      PICOMP(true,Don,Forward, IntrF1+((position-k+3) % 3));
    }
    // On reste intronique
    LBP[IntronF3TA].Update(Data.sig[DATA::Acc].weight[Signal::ForwardNo]);
    INSERTANDPAYTHESLOPE(IntronF3TA,DATA::IntronF);
    */
    // ----------------------------------------------------------------
    // ----------------- Introns de phase -k reverse ------------------
    // ----------------------------------------------------------------
    for (k = 0; k<3; k++) {
      maxi = NINFINITY;
      
      // On quitte un Intr ou un Init reverse
      // TODO: no spliceable stop
      PICOMP(true,Don,Reverse, IntrR1+((Data_Len-position-1-k)%3));
      PICOMP(true,Don,Reverse, InitR1+((Data_Len-position-1-k)%3));

      // On reste intronique
      LBP[IntronR1+k].Update(Data.sig[DATA::Acc].weight[Signal::ReverseNo]);

      INSERTANDPAYTHESLOPE(IntronR1+k,DATA::IntronR);
    }
    // ----------------------------------------------------------------
    // ----------------- Introns reverse speciaux ---------------------
    // ----------------------------------------------------------------
    //
    // --- Intron Phase 1 after a G (AT)
    //
    /* TODO spliceable STOPs
       maxi = NINFINITY;
    k = 2;
    if (StopStop & DNASeq::StopAfterG) {
      // - on quitte un Intr ou un Term
      PICOMP(true,Acc,Reverse, IntrR1+((Data_Len-position-k) % 3));
      PICOMP(true,Acc,Reverse, TermR1+((Data_Len-position-k) % 3));
    }
    // On reste intronique
    LBP[IntronR3G].Update(Data.sig[DATA::Acc].weight[Signal::ReverseNo]);
    INSERTANDPAYTHESLOPE(IntronR3G,DATA::IntronR);

    //
    // --- Intron Phase 1 after an A (GT|AT)
    //
    maxi = NINFINITY;
    k = 2;
    if (StopStop & DNASeq::StopAfterA) {
      // - on quitte un Intr ou un Term
      PICOMP(true,Acc,Reverse, IntrR1+((Data_Len-position-k) % 3));
      PICOMP(true,Acc,Reverse, TermR1+((Data_Len-position-k) % 3));
    }
    // On reste intronique
    LBP[IntronR3A].Update(Data.sig[DATA::Acc].weight[Signal::ReverseNo]);
    INSERTANDPAYTHESLOPE(IntronR3A,DATA::IntronR);

    //
    // --- Intron Phase 2 after an AG, AA ou GA (T)
    //
    maxi = NINFINITY;
    k = 1;
    if (StopStop & DNASeq::StopAfterAG) {
      // - on quitte un Intr ou un Term
      PICOMP(true,Acc,Reverse, IntrR1+((Data_Len-position-k) % 3));
      PICOMP(true,Acc,Reverse, TermR1+((Data_Len-position-k) % 3));
    }
    // On reste intronique
    LBP[IntronR2AG].Update(Data.sig[DATA::Acc].weight[Signal::ReverseNo]);
    INSERTANDPAYTHESLOPE(IntronR2AG,DATA::IntronR);
    */
}
