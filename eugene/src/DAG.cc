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

//#define DEBUGME 1

extern Parameters PAR;

// ----------------------------------------------------------------
// Default constructor.
// ----------------------------------------------------------------
DAG :: DAG ()
{
  StartPosition=0;
  EndPosition=0;
  EvidenceName[0]=0;
  pred = NULL;
  NormalizingPath= 0.0;
  TheSeq = NULL;

  for (int i = 0; i < NbTracks; i++) 
    LBP[i].Path.InitState(i,StartPosition);
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
  for (int i = 0; i < NbTracks; i++) LBP[i].Path.InitState(i,StartPosition);
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
   for (int i = 0; i < NbTracks; i++) LBP[i].Path.InitState(i,StartPosition);
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
  
  int j,k = (Forward ? TheSeq->SeqLen+1 :-1);

  // Insert best possible backpoint at the start of the algo (the
  // insert is not automatically possible, cf cost dist. 
  // on length)
  for (j = 0; j < NbTracks; j++) {
    PrevBP[j] = LBP[j].BestUsable(k,&PBest[j],0);
    LBP[j].ForceNew(j,k,PBest[j],PrevBP[j]);
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
  pred->optimalPath = maxi+NormalizingPath;
  return maxi+NormalizingPath;
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
//  Print active BP stats
// ----------------------------------------------------------------
void DAG :: StatActive() {
  int BPalloc = 0;

  for (int k = 0;k < NbTracks ; k++) 
    BPalloc += LBP[k].NumBPAlloc;
  
   printf("Number of active BP: %d\n",BPalloc);
}
// ----------------------------------------------------------------
//  Print collected BP stats
// ----------------------------------------------------------------
void DAG :: StatGC() {
  int BPcollect = 0;

  for (int k = 0;k < NbTracks ; k++) 
    BPcollect += LBP[k].NumBPCollect;
  
   printf("Number of collected BP: %d\n",BPcollect);
}
// ----------------------------------------------------------------
//  MarkandSweep garbage collector
// ----------------------------------------------------------------
void DAG :: MarkAndSweep(int pos){
  
  int k;
  int Horizon = pos;

  printf("GC started, ");
  StatActive();

  // We first compute the maximum horizon that will be GC'd. This is
  // a purely approximative approach. There is a priori non guarantee
  // that everything will be GC'd with this but most things should.
  // For each track, the penalty distribution may have a different
  // length and everything inside this length is a source pointeur for
  // tracing. The horizon is therefore the Minimum over all Tracks, of
  // the maximum of 0 and (currentpos - PenD.MaxLen - k*GCLATENCY). k
  // = 2 should be nice.

  
  for (k = 0; k < NbTracks ; k++) 
    Horizon = Min(Horizon,Max(0,pos - LBP[k].PenD.MaxLen - k*GCLATENCY));

  for (k = 0; k < NbTracks ; k++)
    LBP[k].ClearMark(Horizon);

  for (k = 0; k < NbTracks ; k++) 
    LBP[k].Mark(Horizon);
  
  for (k = 0; k < NbTracks ; k++) 
    LBP[k].Sweep(Horizon);

  printf("GC finished, ");
  StatGC();
}

// ----------------------------------------------------------------
//  Normalize to avoid rounding errors
// ----------------------------------------------------------------
inline void DAG::Normalize()
{
  double BestU,maxi = -NINFINITY;

  for (int k = 0 ; k < NbTracks; k++) {
    BestU = LBP[k].Optimal;
    if ((BestU > NINFINITY) && (BestU < maxi)) 
      maxi = BestU;
  }
  
  for (int k = 0 ; k < NbTracks; k++) 
    LBP[k].Update(-maxi);
  NormalizingPath += maxi;
}

// ----------------------------------------------------------------
//  Account for length penalty distributions
// ----------------------------------------------------------------
inline void DAG::ApplyLengthPenalty(int position, DATA Data, int NoContentsUpdate)
{
  int Data_Len = TheSeq->SeqLen;
    
  if ((position < Data_Len) && (!NoContentsUpdate))
    for (int k=0; k < NbTracks; k++)
      LBP[k].PayTheSlope();
}

// ----------------------------------------------------------------
//  Account for contents score
// ----------------------------------------------------------------
void DAG::ApplyScore(int position, DATA Data, int NoContentsUpdate)
{
  int Data_Len = TheSeq->SeqLen;
  
  if ((position < Data_Len) && (!NoContentsUpdate)){

    for (int k=0; k< NbTracks; k++)
      LBP[k].Update(Data.contents[SensorContents[k]]);
  }
}
// ----------------------------------------------------------------
// A first set of macros
// ----------------------------------------------------------------
#define ISPOSSIBLE(X,Y) (!(isinf(Data.sig[DATA::X].weight[Signal::Y])))

#define INEED(K) if (!PrevBP[K])					\
    (PrevBP[Strand ?  ReverseIt[K] : K] =				\
     LBP[Strand ? ReverseIt[K] : K].BestUsable(position, &PBest[Strand ? ReverseIt[K] : K]))

// ----------------------------------------------------------------
// Calcul des meilleures opening edges
// ----------------------------------------------------------------
inline void DAG::ComputeRequired(enum Signal::Edge Strand, DATA Data, int position)
{
  int Data_Len = TheSeq->SeqLen;
  int PhaseF = (Strand ? ((Data_Len-position+3) % 3) : (position % 3));
  int PhaseR = (Strand ? (position % 3) : ((Data_Len-position+3) % 3));

  for (int k = 0; k < NbTracks; k++) {
    PrevBP[k] = NULL;
    PBest[k] = NINFINITY;
  }
  
  // ---------------------------
  // ExonsR spliced (from Intron)
  // ---------------------------
  if (ISPOSSIBLE(Don,Reverse-Strand)) {
    INEED(IntronR1); 
    INEED(IntronR2); 
    INEED(IntronR3);
    INEED(IntronR2AG);
    INEED(IntronR3G);
    INEED(IntronR3A);
  }
  
  // ---------------------------
  // ExonsR (from UTR3R)
  // ---------------------------
  if (ISPOSSIBLE(Stop,Reverse-Strand)) 
    INEED(UTR3R);
    
  // ---------------------------
  // ExonsR (from frameshift)
  // ---------------------------
  if (ISPOSSIBLE(Ins,Reverse-Strand) || ISPOSSIBLE(Del,Reverse-Strand)) {
    INEED(InitR1); INEED(InitR2); INEED(InitR3);
    INEED(IntrR1); INEED(IntrR2); INEED(IntrR3);
    INEED(TermR1); INEED(TermR2); INEED(TermR3); 
    INEED(SnglR1); INEED(SnglR2); INEED(SnglR3); 
  }
  
  // ---------------------------
  // IntronR (from ExonR)
  // ---------------------------
  if (ISPOSSIBLE(Acc,Reverse-Strand)) {
    INEED(IntrR1); INEED(IntrR2); INEED(IntrR3);
    INEED(TermR1); INEED(TermR2); INEED(TermR3); 
  }

  // ---------------------------
  // UTR5R (from Init/SnglR) - check
  // ---------------------------
  if (ISPOSSIBLE(Start,Reverse-Strand)) {
    INEED(InitR1+PhaseR);
    INEED(SnglR1+PhaseR); 
  }

  // ---------------------------
  // UTR5R (from IntronU5R)
  // ---------------------------
  if (ISPOSSIBLE(Don,Reverse-Strand))
    INEED(IntronU5R);

  // ---------------------------
  // IntronU5R (from UTR5R)
  // ---------------------------
  if (ISPOSSIBLE(Acc,Reverse-Strand))
    INEED(UTR5R);

  // ---------------------------
  // UTR3R (from InterGen)
  // ---------------------------
  if (ISPOSSIBLE(tStop,Reverse-Strand))
    INEED(InterGen);

  // ---------------------------
  // UTR3R (from IntronU3R)
  // ---------------------------
  if (ISPOSSIBLE(Don,Reverse-Strand))
    INEED(IntronU3R);

  // ---------------------------
  // IntronU3R (from UTR3R)
  // ---------------------------
  if (ISPOSSIBLE(Acc,Reverse-Strand))
    INEED(UTR3R);

  // ---------------------------
  // InterG (from UTR5R)
  // ---------------------------
  if (ISPOSSIBLE(tStart,Reverse-Strand))
    INEED(UTR5R);

  // ---------------------------
  // InterG (from UTR3F)
  // ---------------------------
  if (ISPOSSIBLE(tStop,Forward+Strand))
    INEED(UTR3F);

  // ---------------------------
  // IntronU5F (from UTR5F)
  // ---------------------------
  if (ISPOSSIBLE(Don,Forward+Strand))
    INEED(UTR5F);

  // ---------------------------
  // UTR5F (from InterGen)
  // ---------------------------
  if (ISPOSSIBLE(tStart,Forward+Strand))
    INEED(InterGen);

  // ---------------------------
  // UTR5F (from IntronU5F)
  // ---------------------------
  if (ISPOSSIBLE(Acc,Forward+Strand))
    INEED(IntronU5F);

  // ---------------------------
  // IntronU3R (from UTR3F)
  // ---------------------------
  if (ISPOSSIBLE(Don,Forward+Strand))
    INEED(UTR3F);

  // ---------------------------
  // UTR3F (from IntronU3F)
  // ---------------------------
  if (ISPOSSIBLE(Acc,Forward+Strand))
    INEED(IntronU3F);

  // ---------------------------
  // UTR3F (from ExonF)
  // ---------------------------
  if (ISPOSSIBLE(Stop,Forward+Strand)) {
    INEED(TermF1+PhaseF);
    INEED(SnglF1+PhaseF); 
  }

  // ---------------------------
  // ExonF (from UTR5F)
  // ---------------------------
  if (ISPOSSIBLE(Start,Forward+Strand))
    INEED(UTR5F);

  // ---------------------------
  // ExonF (from IntronF)
  // ---------------------------
  if (ISPOSSIBLE(Acc,Forward+Strand)) {
    INEED(IntronF2T);
    INEED(IntronF3TG);
    INEED(IntronF3TA);      
    INEED(IntronF1); 
    INEED(IntronF2); 
    INEED(IntronF3);
  }

  // ---------------------------
  // ExonF (from frameshift)
  // ---------------------------
  if (ISPOSSIBLE(Ins,Forward+Strand) || 
      ISPOSSIBLE(Del,Forward+Strand)) {
    INEED(InitF1); INEED(InitF2); INEED(InitF3);
    INEED(IntrF1); INEED(IntrF2); INEED(IntrF3);
    INEED(TermF1); INEED(TermF2); INEED(TermF3); 
    INEED(SnglF1); INEED(SnglF2); INEED(SnglF3); 
  }

  // ---------------------------
  // IntronF (from ExonF)
  // ---------------------------
  if (ISPOSSIBLE(Don,Forward+Strand)) {
    INEED(InitF1); INEED(InitF2); INEED(InitF3);
    INEED(IntrF1); INEED(IntrF2); INEED(IntrF3);
  }

#ifdef DEBUGME
    for (k=0; k<(NbTracks-1)/2;k++)
      if (PrevBP[ReverseTracks[k]]) 
	printf("Cost %f Track %d a From %d\n",
	       PBest[ReverseTracks[k]],k,ReverseIt[PrevBP[ReverseTracks[k]]->State]);

      if (PrevBP[UnorientedTracks[0]]) 
	printf("Cost %f Track %d From %d\n",
	       PBest[UnorientedTracks[0]],k,ReverseIt[PrevBP[UnorientedTracks[0]]->State]);

    for (k=0; k<(NbTracks-1)/2;k++)
      if (PrevBP[ForwardTracks[k]]) 
	printf("Cost %f Track %d b From %d\n",
	       PBest[ForwardTracks[k]],k,ReverseIt[PrevBP[ForwardTracks[k]]->State]);

    printf("---------- pos %d/%d, best %d norm %f ----------\n",position,Data_Len-position,ReverseIt[(int)best],NormalizingPath);
#endif
}
// ----------------------------------------------------------------
// A second set of macros
// ----------------------------------------------------------------
#define PICOMPEN(C,S,B,O,P) if ((C) && ISPOSSIBLE(S,B)) {               \
    BestU = PBest[Strand ? ReverseIt[O] :O]+Data.sig[DATA::S].weight[Signal::B]+(P); \
    if (isnan(maxi) || (BestU > maxi)) {maxi = BestU; best = (Strand ? ReverseIt[O] :O);}}

#define PICOMP(C,S,B,O) PICOMPEN(C,S,B,O,0.0)

#define INSERT(P)					\
  LBP[Strand ? ReverseIt[P] :P].InsertNew(best, position, maxi,PrevBP[best]);

// ----------------------------------------------------------------
// Perform DP itself: one recursive level through possible signals
// ----------------------------------------------------------------
inline void DAG::ComputeSigShifts(enum Signal::Edge Strand, DATA Data, int position)
{
  int k;
  int Data_Len = TheSeq->SeqLen;
  int PhaseF = (Strand ? ((Data_Len-position+3) % 3) : (position % 3));
  int PhaseR = (Strand ? (position % 3) : ((Data_Len-position+3) % 3));
  double BestU, maxi = -NINFINITY;
  signed   char best = 'Z';

  // Get information on possible spliced stops
  int StopStop = TheSeq->IsStopStop(position);
  int StartStop = TheSeq->IsStartStop(position);

  // ----------------------------------------------------------------
  // ------------------ Inits en forward ----------------------------
  // ----------------------------------------------------------------
  for (k = 0; k < 3; k++) {
    maxi = NINFINITY; best = -1;     
    
    // On commence a coder (Start),si en phase ca vient d'une UTR 5' forward
    PICOMP((PhaseF == k),Start,Forward+Strand,UTR5F);
    // Il y a une insertion (frameshift). Saut de position de nucléotide ignoré.
    PICOMP(true,Ins,Forward+Strand,InitF1+(k+1)%3);
    // Il y a une deletion (frameshift)
    PICOMP(true,Del,Forward+Strand,InitF1+(k+2)%3);
    
    // On va tout droit.
      // S'il y  a un STOP en phase on ne peut continuer
    if (PhaseF == k)
      LBP[InitF1+k].Update(Data.sig[DATA::Stop].weight[Signal::ForwardNo+Strand]);
    LBP[InitF1+k].Update(Data.sig[DATA::Don].weight[Signal::ForwardNo+Strand]);
    
    INSERT(InitF1+k);
  }
  // ----------------------------------------------------------------
  // ------------------------- Inits en reverse ---------------------
  // ----------------------------------------------------------------
  for (k = 3; k<6; k++) {
    maxi = NINFINITY; best = -1;
    
    // - on recommence a coder (Donneur) Ca vient d'un intron (no spliceable stop)
    PICOMPEN(true,Don,Reverse-Strand,IntronR1+((Data_Len-position-k) % 3),
	     ((PhaseR == k-3) ? Data.sig[DATA::Stop].weight[Signal::ReverseNo-Strand] : 0.0));
    
    // Not AfterG
    PICOMPEN((Data_Len-position-k) % 3 == 2,Don,Reverse-Strand,IntronR3G,
	     ((StartStop & DNASeq::isTAr) ? SplicedStopPen : 0.0));
    // Not AfterA
    PICOMPEN((Data_Len-position-k) % 3 == 2,Don,Reverse-Strand,IntronR3A ,
	     ((StartStop & (DNASeq::isTGr | DNASeq::isTAr))  ? SplicedStopPen : 0.0));
    // Not AfterAG
    PICOMPEN((Data_Len-position-k) % 3 == 1,Don,Reverse-Strand,IntronR2AG,
	     ((StartStop & DNASeq::isTr) ? SplicedStopPen : 0.0));
    
    // Il y a une insertion (frameshift)
    PICOMP(true,Ins,Reverse-Strand, InitR1+(k+2)%3);
    // Il y a une deletion (frameshift)
    PICOMP(true,Del,Reverse-Strand, InitR1+(k+1)%3);
    
    // On va tout droit.
    // S'il y  a un STOP en phase on ne peut continuer
    if ((PhaseR % 3 == k-3)) 
      LBP[InitF1+k].Update(Data.sig[DATA::Stop].weight[Signal::ReverseNo-Strand]);
    LBP[InitF1+k].Update(Data.sig[DATA::Don].weight[Signal::ReverseNo-Strand]);
    
    INSERT(InitF1+k);
  }
  // ----------------------------------------------------------------
  // ------------------ Sngl en forward ----------------------------
  // ----------------------------------------------------------------
  for (k = 0; k < 3; k++) {
    maxi = NINFINITY; best = -1;     
    
    // On commence a coder (Start),si en phase ca vient d'une UTR 5' forward
    PICOMP((PhaseF == k),Start,Forward+Strand,UTR5F);
    // Il y a une insertion (frameshift). Saut de positionléotide ignore.
    PICOMP(true,Ins,Forward+Strand,SnglF1+(k+1)%3);
    // Il y a une deletion (frameshift)
    PICOMP(true,Del,Forward+Strand,SnglF1+(k+2)%3);
    
    // On va tout droit.
    // S'il y  a un STOP en phase on ne peut continuer
    if ((PhaseF == k))
      LBP[SnglF1+k].Update(Data.sig[DATA::Stop].weight[Signal::ForwardNo+Strand]);
    LBP[SnglF1+k].Update(Data.sig[DATA::Don].weight[Signal::ForwardNo+Strand]);
    
    INSERT(SnglF1+k);
  }
  // ----------------------------------------------------------------
  // ------------------------- Sngl en reverse ---------------------
  // ----------------------------------------------------------------
  for (k = 3; k<6; k++) {
    maxi = NINFINITY; best = -1;
    
    // On commence a coder (Stop). Ca vient d'une UTR 3' reverse
    PICOMP((PhaseR == k-3),Stop,Reverse-Strand,UTR3R);
    // Il y a une insertion (frameshift)
    PICOMP(true,Ins,Reverse-Strand, SnglR1+(k+2)%3);
    // Il y a une deletion (frameshift)
    PICOMP(true,Del,Reverse-Strand, SnglR1+(k+1)%3);
    
    // On va tout droit.
    // S'il y  a un STOP en phase on ne peut continuer
    if (PhaseR == k-3)
      LBP[SnglF1+k].Update(Data.sig[DATA::Stop].weight[Signal::ReverseNo-Strand]);
    LBP[SnglF1+k].Update(Data.sig[DATA::Don].weight[Signal::ReverseNo-Strand]);
    
    INSERT(SnglF1+k);
  }
  // ----------------------------------------------------------------
  // ------------------ Intrs en forward ----------------------------
  // ----------------------------------------------------------------
  for (k = 0; k < 3; k++) {
    maxi = NINFINITY; best = -1;     
    
    // On recommence a coder (Accepteur). Ca vient d'un intron (no spliceable stop)
    PICOMP(true,Acc,Forward+Strand,IntronF1+((position-k+3) % 3));
    
    // Not AfterT
    PICOMPEN(((position-k+3) % 3) == 1,Acc,Forward+Strand,IntronF2T ,
	     ((StopStop & (DNASeq::isGAf | DNASeq::isARf))  ? SplicedStopPen : 0.0));
    // Not AfterTG
    PICOMPEN(((position-k+3) % 3) == 2,Acc,Forward+Strand,IntronF3TG,
	     ((StopStop & DNASeq::isAf) ? SplicedStopPen : 0.0));
    // Not AfterTA
    PICOMPEN(((position-k+3) % 3) == 2,Acc,Forward+Strand,IntronF3TA,
	     ((StopStop & (DNASeq::isGf | DNASeq::isAf)) ? SplicedStopPen : 0.0));
    
    // Il y a une insertion (frameshift). Saut de position de nucléotide ignoré.
    PICOMP(true,Ins,Forward+Strand,IntrF1+(k+1)%3);
    // Il y a une deletion (frameshift)
    PICOMP(true,Del,Forward+Strand,IntrF1+(k+2)%3);
    
    // On va tout droit.
    // S'il y  a un STOP en phase on ne peut continuer
    if (PhaseF == k)
      LBP[IntrF1+k].Update(Data.sig[DATA::Stop].weight[Signal::ForwardNo+Strand]);
    LBP[IntrF1+k].Update(Data.sig[DATA::Don].weight[Signal::ForwardNo+Strand]);
    
    INSERT(IntrF1+k);
  }
  // ----------------------------------------------------------------
  // ------------------------- Intrs en reverse ---------------------
  // ----------------------------------------------------------------
  for (k = 3; k<6; k++) {
    maxi = NINFINITY; best = -1;
    
    // - on recommence a coder (Donneur) Ca vient d'un intron (no spliceable stop)
    PICOMPEN(true,Don,Reverse-Strand,IntronR1+((Data_Len-position-k) % 3),((PhaseR == k-3) ? Data.sig[DATA::Stop].weight[Signal::ReverseNo-Strand] : 0.0));
    
    // Not AfterG
    PICOMPEN(((Data_Len-position-k) % 3) == 2,Don,Reverse-Strand,IntronR3G ,
	     ((StartStop & DNASeq::isTAr)  ? SplicedStopPen : 0.0));
    // Not AfterA
    PICOMPEN(((Data_Len-position-k) % 3) == 2,Don,Reverse-Strand,IntronR3A ,
	     ((StartStop & (DNASeq::isTAr | DNASeq::isTGr))  ? SplicedStopPen : 0.0));
    // Not AfterAG
    PICOMPEN(((Data_Len-position-k) % 3) == 1,Don,Reverse-Strand,IntronR2AG,
	     ((StartStop & DNASeq::isTr) ? SplicedStopPen : 0.0));
    
    // Il y a une insertion (frameshift)
    PICOMP(true,Ins,Reverse-Strand, IntrR1+(k+2)%3);
    // Il y a une deletion (frameshift)
    PICOMP(true,Del,Reverse-Strand, IntrR1+(k+1)%3);
    
    // On va tout droit.
    // S'il y  a un STOP en phase on ne peut continuer
    if (PhaseR == k-3) 
      LBP[IntrF1+k].Update(Data.sig[DATA::Stop].weight[Signal::ReverseNo-Strand]);
    LBP[IntrF1+k].Update(Data.sig[DATA::Don].weight[Signal::ReverseNo-Strand]);
    
    INSERT(IntrF1+k);
  }
  // ----------------------------------------------------------------
  // ------------------ Terms en forward ----------------------------
  // ----------------------------------------------------------------
  for (k = 0; k < 3; k++) {
    maxi = NINFINITY; best = -1;     
    
    // On recommence a coder (Accepteur). Ca vient d'un intron
    PICOMP(true,Acc,Forward+Strand,IntronF1+((position-k+3) % 3));
    
    // Not AfterT
    PICOMPEN(((position-k+3) % 3) == 1,Acc,Forward+Strand,IntronF2T ,
	     ((StopStop & (DNASeq::isGAf | DNASeq::isARf))  ? SplicedStopPen :0.0));
    // Not AfterTG
    PICOMPEN(((position-k+3) % 3) == 2,Acc,Forward+Strand,IntronF3TG,
	     ((StopStop & DNASeq::isAf) ? SplicedStopPen : 0.0));
    // Not AfterTA
    PICOMPEN(((position-k+3) % 3) == 2,Acc,Forward+Strand,IntronF3TA,
	     ((StopStop & (DNASeq::isAf | DNASeq::isGf)) ? SplicedStopPen : 0.0));
    
    // Il y a une insertion (frameshift). Saut de positionléotide ignore.
    PICOMP(true,Ins,Forward+Strand,TermF1+(k+1)%3);
    // Il y a une deletion (frameshift)
    PICOMP(true,Del,Forward+Strand,TermF1+(k+2)%3);
    
    // On va tout droit.
    // S'il y  a un STOP en phase on ne peut continuer
    if (PhaseF == k)
      LBP[TermF1+k].Update(Data.sig[DATA::Stop].weight[Signal::ForwardNo+Strand]);
    LBP[TermF1+k].Update(Data.sig[DATA::Don].weight[Signal::ForwardNo+Strand]);
    
    INSERT(TermF1+k);
  }
  // ----------------------------------------------------------------
  // ------------------------- Terms en reverse ---------------------
  // ----------------------------------------------------------------
  for (k = 3; k<6; k++) {
    maxi = NINFINITY; best = -1;
    
    // On commence a coder (Stop). Ca vient d'une UTR 3' reverse
    PICOMP((PhaseR == k-3),Stop,Reverse-Strand,UTR3R);
    // Il y a une insertion (frameshift)
    PICOMP(true,Ins,Reverse-Strand, TermR1+(k+2)%3);
    // Il y a une deletion (frameshift)
    PICOMP(true,Del,Reverse-Strand, TermR1+(k+1)%3);
    
    // On va tout droit.
    // S'il y  a un STOP en phase on ne peut continuer
    if ((PhaseR == k-3)) 
      LBP[TermF1+k].Update(Data.sig[DATA::Stop].weight[Signal::ReverseNo-Strand]);
    LBP[TermF1+k].Update(Data.sig[DATA::Don].weight[Signal::ReverseNo-Strand]);
    
    INSERT(TermF1+k);
  }
  // ----------------------------------------------------------------
  // ------------------------ Intergenique --------------------------
  // ----------------------------------------------------------------
  // Ca peut venir d'une fin de 3' direct ou de 5' reverse
  maxi = NINFINITY; best = -1;
  
  // From 5' reverse
  PICOMP(true,tStart,Reverse-Strand, UTR5R);
  // From 3' direct
  PICOMP(true,tStop,Forward+Strand, UTR3F);
  
  // On reste intergenique
  // et les transstartNO/TransstopNO ???
  
  INSERT(InterGen);
  // ----------------------------------------------------------------
  // ---------------------- UTR 5' direct ---------------------------
  // ----------------------------------------------------------------
  maxi = NINFINITY; best = -1;
  
  // On vient de l'intergenique. 
  PICOMP(true,tStart,Forward+Strand, InterGen);
  // On peut venir aussi d'un intron d'UTR.
  PICOMP(true,Acc,Forward+Strand, IntronU5F);
  
  // On reste 5' direct. On ne prend pas le Start eventuel.
  LBP[UTR5F].Update(Data.sig[DATA::Start].weight[Signal::ForwardNo+Strand]);
  
  INSERT(UTR5F);
  // ----------------------------------------------------------------
  // ------------------- Introns d'UTR5F ----------------------------
  // ----------------------------------------------------------------
  maxi = NINFINITY; best = -1;
  
  // on quitte l'UTR 
  PICOMP(true,Don,Forward+Strand, UTR5F);
  
  // On reste intronique
  LBP[IntronU5F].Update(Data.sig[DATA::Acc].weight[Signal::ForwardNo+Strand]);
  
  INSERT(IntronU5F);
  // ----------------------------------------------------------------
  // ------------------- Introns d'UTR3F ----------------------------
  // ----------------------------------------------------------------
  maxi = NINFINITY; best = -1;
  
  // on quitte l'UTR 
  PICOMP(true,Don,Forward+Strand, UTR3F);
  
  // On reste intronique
  LBP[IntronU3F].Update(Data.sig[DATA::Acc].weight[Signal::ForwardNo+Strand]);
  
  INSERT(IntronU3F);
  // ----------------------------------------------------------------
  // ---------------------- UTR 3' direct ---------------------------
  // ----------------------------------------------------------------
  maxi = NINFINITY; best = -1;
  
  // Ca vient d'un Term ou d'un Sngl direct + STOP
  PICOMP(true,Stop,Forward+Strand, TermF1+PhaseF);
  PICOMP(true,Stop,Forward+Strand, SnglF1+PhaseF);
  // On peut venir aussi d'un intron d'UTR.
  PICOMP(true,Acc,Forward+Strand, IntronU3F);
  
  // On reste 3' direct
  
  INSERT(UTR3F);
  // ----------------------------------------------------------------
  // ----------------------- UTR 5'reverse --------------------------
  // ----------------------------------------------------------------
  maxi = NINFINITY; best = -1;
  
  // Ca vient d'un Sngl ou Init reverse + START
  PICOMP(true,Start,Reverse-Strand, InitR1+PhaseR);
  PICOMP(true,Start,Reverse-Strand, SnglR1+PhaseR);
  // On peut venir aussi d'un intron d'UTR.
  PICOMP(true,Don,Reverse-Strand, IntronU5R);
  
  // On reste 5' reverse
  LBP[UTR5R].Update(Data.sig[DATA::Start].weight[Signal::ReverseNo-Strand]);
  
  INSERT(UTR5R);
  // ----------------------------------------------------------------
  // ------------------- Introns d'UTR5R ----------------------------
  // ----------------------------------------------------------------
  maxi = NINFINITY; best = -1;
  
  // On quitte une UTR5R
  PICOMP(true,Acc,Reverse-Strand, UTR5R);
  
  // On reste intronique
  LBP[IntronU5R].Update(Data.sig[DATA::Acc].weight[Signal::ReverseNo-Strand]);
  
  INSERT(IntronU5R);
  // ----------------------------------------------------------------
  // ------------------- Introns d'UTR3R ----------------------------
  // ----------------------------------------------------------------
  maxi = NINFINITY; best = -1;
  
  // On quitte une UTR5R
  PICOMP(true,Acc,Reverse-Strand, UTR3R);
  
  // On reste intronique
  LBP[IntronU3R].Update(Data.sig[DATA::Acc].weight[Signal::ReverseNo-Strand]);
  
  INSERT(IntronU3R);
  // ----------------------------------------------------------------
  // ----------------------- UTR 3'reverse --------------------------
  // ----------------------------------------------------------------
  maxi = NINFINITY; best = -1;
  
  // On demarre depuis l'intergenique
  PICOMP(true,tStop,Reverse-Strand, InterGen);
  // On peut venir aussi d'un intron d'UTR.
  PICOMP(true,Don,Reverse-Strand, IntronU3R);
  
  // On reste 3' reverse
  
  INSERT(UTR3R);
  // ----------------------------------------------------------------
  // ---------------- Introns de phase k forward --------------------
  // ----------------------------------------------------------------
  for (k = 0; k<3; k++) {
    maxi = NINFINITY; best = -1;
    
    // - on quitte un Init ou un Intr
    // no spliceable stop: 
    if (!(((StartStop & DNASeq::isTf) && k == 1) ||
	  ((StartStop & (DNASeq::isTGf|DNASeq::isTAf)) && k == 2))) {
      PICOMPEN(true,Don,Forward+Strand, InitF1+((position-k+3) % 3),
	       (k == 0 ? Data.sig[DATA::Stop].weight[Signal::ForwardNo+Strand] : 0.0));
      PICOMPEN(true,Don,Forward+Strand, IntrF1+((position-k+3) % 3),
	       (k == 0 ? Data.sig[DATA::Stop].weight[Signal::ForwardNo+Strand] : 0.0));
    }
    
    // On reste intronique
    LBP[IntronF1+k].Update(Data.sig[DATA::Acc].weight[Signal::ForwardNo+Strand]);
    
    INSERT(IntronF1+k);
  }
  // ----------------------------------------------------------------
  // ----------------- Introns forward speciaux ---------------------
  // ----------------------------------------------------------------
  //
  // --- Intron Phase 1 after a T (GA|AA|AG)
  //
  maxi = NINFINITY; best = -1;
  k = 1;
  if (StartStop & DNASeq::isTf) {
    // - on quitte un Init ou un Intr
    PICOMP(true,Don,Forward+Strand, InitF1+((position-k+3) % 3));
    PICOMP(true,Don,Forward+Strand, IntrF1+((position-k+3) % 3));
  }
  // On reste intronique
  LBP[IntronF2T].Update(Data.sig[DATA::Acc].weight[Signal::ForwardNo+Strand]);
  INSERT(IntronF2T);
  
  //
  // --- Intron Phase 2 after an TG|A
  //
  maxi = NINFINITY; best = -1;
  k = 2;
  if (StartStop & DNASeq::isTGf) {
    // - on quitte un Init ou un Intr
    PICOMP(true,Don,Forward+Strand, InitF1+((position-k+3) % 3));
    PICOMP(true,Don,Forward+Strand, IntrF1+((position-k+3) % 3));
  }
  // On reste intronique
  LBP[IntronF3TG].Update(Data.sig[DATA::Acc].weight[Signal::ForwardNo+Strand]);
  INSERT(IntronF3TG);
  
  //
  // --- Intron Phase 2 after a TA(A|G)
  //
  maxi = NINFINITY; best = -1;
  k = 2;
  if (StartStop & DNASeq::isTAf) {
    // - on quitte un Init ou un Intr
    PICOMP(true,Don,Forward+Strand, InitF1+((position-k+3) % 3));
    PICOMP(true,Don,Forward+Strand, IntrF1+((position-k+3) % 3));
  }
  // On reste intronique
  LBP[IntronF3TA].Update(Data.sig[DATA::Acc].weight[Signal::ForwardNo+Strand]);
  INSERT(IntronF3TA);
  // ----------------------------------------------------------------
  // ----------------- Introns de phase -k reverse ------------------
  // ----------------------------------------------------------------
  for (k = 0; k<3; k++) {
    maxi = NINFINITY; best = -1;
    
    // On quitte un Intr ou un Term reverse
    // no spliceable stop: 
    if (!(((StopStop & (DNASeq::isGr|DNASeq::isAr)) && k == 2) ||
	  ((StopStop & (DNASeq::isGAr | DNASeq::isARr)) && k == 1)))  {
      PICOMP(true,Acc,Reverse-Strand, IntrR1+((Data_Len-position-k) % 3));
      PICOMP(true,Acc,Reverse-Strand, TermR1+((Data_Len-position-k) % 3));
    }
    
    // On reste intronique
    LBP[IntronR1+k].Update(Data.sig[DATA::Acc].weight[Signal::ReverseNo-Strand]);
    
    INSERT(IntronR1+k);
  }
  // ----------------------------------------------------------------
  // ----------------- Introns reverse speciaux ---------------------
  // ----------------------------------------------------------------
  //
  // --- Intron Phase 1 after a G (AT)
  //
  maxi = NINFINITY; best = -1;
  k = 2;
  if (StopStop & DNASeq::isGr) {
    // - on quitte un Intr ou un Term
    PICOMP(true,Acc,Reverse-Strand, IntrR1+((Data_Len-position-k) % 3));
    PICOMP(true,Acc,Reverse-Strand, TermR1+((Data_Len-position-k) % 3));
  }
  // On reste intronique
  LBP[IntronR3G].Update(Data.sig[DATA::Acc].weight[Signal::ReverseNo-Strand]);
  INSERT(IntronR3G);
  
  //
  // --- Intron Phase 1 after an A (GT|AT)
  //
  maxi = NINFINITY; best = -1;
  k = 2;
  if (StopStop & DNASeq::isAr) {
    // - on quitte un Intr ou un Term
    PICOMP(true,Acc,Reverse-Strand, IntrR1+((Data_Len-position-k) % 3));
    PICOMP(true,Acc,Reverse-Strand, TermR1+((Data_Len-position-k) % 3));
  }
  // On reste intronique
  LBP[IntronR3A].Update(Data.sig[DATA::Acc].weight[Signal::ReverseNo-Strand]);
  INSERT(IntronR3A);
  
  //
  // --- Intron Phase 2 after an AG, AA ou GA (T)
  //
  maxi = NINFINITY; best = -1;
  k = 1;
  if (StopStop & (DNASeq::isGAr | DNASeq::isARr)) {
    // - on quitte un Intr ou un Term
    PICOMP(true,Acc,Reverse-Strand, IntrR1+((Data_Len-position-k) % 3));
    PICOMP(true,Acc,Reverse-Strand, TermR1+((Data_Len-position-k) % 3));
  }
  // On reste intronique
  LBP[IntronR2AG].Update(Data.sig[DATA::Acc].weight[Signal::ReverseNo-Strand]);
  INSERT(IntronR2AG); 
}

// ----------------------------------------------------------------
//  Shortest Path Algorithm with length constraint
// ----------------------------------------------------------------
void DAG :: ShortestPathAlgoForward (int position, DATA Data)
{
  // Avoid rounding errors on long sequences
  Normalize();
  // Precompute required track values
  ComputeRequired(Signal::Forward,Data,position);
  //Recursive computation
  ComputeSigShifts(Signal::Forward, Data, position);
  // Account for Score and length penalty 
  ApplyScore(position, Data, 0);
  ApplyLengthPenalty(position, Data, 0);
}

// ----------------------------------------------------------------
//  Shortest Path Algorithm with length constraint (reverse)
// ----------------------------------------------------------------
void DAG :: ShortestPathAlgoBackward (int position, DATA Data, int NoContentsUpdate)
{
  ApplyScore(position, Data, NoContentsUpdate);
  ApplyLengthPenalty(position, Data, NoContentsUpdate);
  
  // Avoid rounding errors on long sequences
  Normalize();
  // Precompute required track values
  ComputeRequired(Signal::Reverse,Data,position);
  //Recursive computation
  ComputeSigShifts(Signal::Reverse, Data, position);
}
