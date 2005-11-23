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

int ReverseTrack(int track)
{
  if (track < 0) return track;
  if (track == InterGen) return InterGen;
  
  // look in forward
  for (int i = 0; i< (NbTracks-1)/2; i++)
    if (ForwardTracks[i] == track) return ReverseTracks[i];

  // look in forward
  for (int i = 0; i< (NbTracks-1)/2; i++)
    if (ReverseTracks[i] == track) return ForwardTracks[i];
  return track;
}

// ----------------------------------------------------------------
// Default constructor.
// ----------------------------------------------------------------
DAG :: DAG ()
{
  StartPosition=0;
    EndPosition=0;
  EvidenceName[0]=0;
  pred = NULL;
  for (int i = 0; i < NbTracks; i++) LBP[i].Path.InitState(i,StartPosition);
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

// My ugly macros

#define ISPOSSIBLE(X,Y) (!(isinf(Data.sig[DATA::X].weight[Signal::Y])))
  //  ((1 && printf("IsPossible(%d,%d)=%d\n",(DATA::X),(Signal::Y),(!(isinf(Data.sig[DATA::X].weight[Signal::Y]))))) ? \
//   (!(isinf(Data.sig[DATA::X].weight[Signal::Y]))) : 0)
#define INEED(K) if (!PrevBP[K]) PrevBP[K] = LBP[K].BestUsable(position, &PBest[K])
  //        printf("INeed BestU[%d]=%f\n",(K),PBest[K])
#define PICOMP(C,S,B,O) if ((C) && ISPOSSIBLE(S,B)) {           \
        BestU = PBest[O]+Data.sig[DATA::S].weight[Signal::B];\
        if (isnan(maxi) || (BestU > maxi)) {maxi = BestU; best = O;}}

  //  	printf("PiComp NewBestU[%d]=%f, Signal[%d][%d]=%f Line %d\n",ReverseTrack(O),BestU+NormalizingPath,DATA::S,1-Signal::B,Data.sig[DATA::S].weight[Signal::B],__LINE__); \


#define PICOMPR(C,S,B,O) if ((C) && ISPOSSIBLE(S,B)) {           \
        BestU = PBest[O]+Data.sig[DATA::S].weight[Signal::B];\
        if (isnan(maxi) || (BestU > maxi)) {maxi = BestU; best = O;}}

  //	printf("PiComp NewBestU[%d]=%f, Signal[%d][%d]=%f Line %d\n",(O),BestU+NormalizingPath,DATA::S,Signal::B,Data.sig[DATA::S].weight[Signal::B],__LINE__); \

#define PICOMPEN(C,S,B,O,P) if ((C) && ISPOSSIBLE(S,B)) {               \
        BestU = PBest[O]+Data.sig[DATA::S].weight[Signal::B]+(P);       \
        if (isnan(maxi) || (BestU > maxi)) {maxi = BestU; best = O;}}

#define INSERT(P)					\
  LBP[P].InsertNew(best, position, maxi,PrevBP[best]);
  
  //  printf("Ins Cost %f From %d Track %d Line %d\n",maxi,best,(P),000); \
  //  printf("Ins Cost %f From %d Track %d Line %d\n",maxi,ReverseTrack(best),(ReverseTrack(P)),000); \


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

  // Avoid rounding errors on long sequences
  Normalize();

  // ----------------------------------------------------------------
  // Calcul des meilleures opening edges
  // ----------------------------------------------------------------
    for (k = 0; k < NbTracks; k++) {
      PrevBP[k] = NULL;
      PBest[k] = NINFINITY;
    }

    // Exons F (No splice)
    if (ISPOSSIBLE(Start,Forward))
      INEED(UTR5F);
    
    if (ISPOSSIBLE(Ins,Forward) || ISPOSSIBLE(Del,Forward)) {
      INEED(InitF1); INEED(InitF2); INEED(InitF3);
      INEED(IntrF1); INEED(IntrF2); INEED(IntrF3);
      INEED(TermF1); INEED(TermF2); INEED(TermF3); 
      INEED(SnglF1); INEED(SnglF2); INEED(SnglF3); 
    }

    // Exons F  (splicing, with and without spliced stops)
    if (ISPOSSIBLE(Acc,Forward)) {
      INEED(IntronF2T);
      INEED(IntronF3TG);
      INEED(IntronF3TA);      
      INEED(IntronF1); 
      INEED(IntronF2); 
      INEED(IntronF3);
    }

    // Exons R  (No Splice)
    if (ISPOSSIBLE(Stop,Reverse)) 
      INEED(UTR3R);
    
    if (ISPOSSIBLE(Ins,Reverse) || ISPOSSIBLE(Del,Reverse)) {
      INEED(InitR1); INEED(InitR2); INEED(InitR3);
      INEED(IntrR1); INEED(IntrR2); INEED(IntrR3);
      INEED(TermR1); INEED(TermR2); INEED(TermR3); 
      INEED(SnglR1); INEED(SnglR2); INEED(SnglR3); 
    }
      
    // Exons R (splicing, with and without spliced stops)
    if (ISPOSSIBLE(Don,Reverse)) {
      INEED(IntronR3G);
      INEED(IntronR3A);
      INEED(IntronR2AG);
      INEED(IntronR1); 
      INEED(IntronR2); 
      INEED(IntronR3);
    }

    // Introns F et R (with no spliceable stops)
    if (ISPOSSIBLE(Don,Forward)) {
      INEED(InitF1); INEED(InitF2); INEED(InitF3);
      INEED(IntrF1); INEED(IntrF2); INEED(IntrF3);
    }

    if (ISPOSSIBLE(Acc,Reverse)) {
      INEED(IntrR1); INEED(IntrR2); INEED(IntrR3);
      INEED(TermR1); INEED(TermR2); INEED(TermR3); 
    }

    // Intergenique
    if (ISPOSSIBLE(tStart,Reverse))
      INEED(UTR5R);
    if (ISPOSSIBLE(tStop,Forward))
      INEED(UTR3F);

    // UTR 5' direct
    if (ISPOSSIBLE(tStart,Forward))
      INEED(InterGen);
    if (ISPOSSIBLE(Acc,Forward))
      INEED(IntronU5F);

    // Introns d'UTR5F
    if (ISPOSSIBLE(Don,Forward))
      INEED(UTR5F);

    // Introns d'UTR3F
    if (ISPOSSIBLE(Don,Forward))
      INEED(UTR3F);

    // UTR 3' direct
    if (ISPOSSIBLE(Acc,Forward))
      INEED(IntronU3F);
    if (ISPOSSIBLE(Stop,Forward)) {
      INEED(TermF1+position%3);
      INEED(SnglF1+position%3); 
    }

    // UTR 5'reverse
    if (ISPOSSIBLE(Start,Reverse)) {
      INEED(InitR1+((Data_Len-position) % 3));
      INEED(SnglR1+((Data_Len-position) % 3)); 
    }

    if (ISPOSSIBLE(Don,Reverse))
      INEED(IntronU5R);

    // Introns d'UTR5R
    if (ISPOSSIBLE(Acc,Reverse))
      INEED(UTR5R);

    // Introns d'UTR5R
    if (ISPOSSIBLE(Acc,Reverse))
      INEED(UTR3R);

    // UTR 3' reverse
    if (ISPOSSIBLE(tStop,Reverse))
      INEED(InterGen);

    if (ISPOSSIBLE(Don,Reverse))
      INEED(IntronU3R);

#ifdef DEBUGME
    for (k=0; k<(NbTracks-1)/2;k++)
      if (PrevBP[ReverseTracks[k]]) 
	printf("Cost %f Track %d a From %d\n",
	       PBest[ReverseTracks[k]],k,ReverseTrack(PrevBP[ReverseTracks[k]]->State));

      if (PrevBP[UnorientedTracks[0]]) 
	printf("Cost %f Track %d From %d\n",
	       PBest[UnorientedTracks[0]],k,ReverseTrack(PrevBP[UnorientedTracks[0]]->State));

    for (k=0; k<(NbTracks-1)/2;k++)
      if (PrevBP[ForwardTracks[k]]) 
	printf("Cost %f Track %d b From %d\n",
	       PBest[ForwardTracks[k]],k,ReverseTrack(PrevBP[ForwardTracks[k]]->State));

    printf("---------- pos %d/%d, best %d norm %f ----------\n",position,Data_Len-position,ReverseTrack((int)best),NormalizingPath);
#endif
      
    // ----------------------------------------------------------------
    // ------------------ Inits en forward ----------------------------
    // ----------------------------------------------------------------
    for (k = 0; k < 3; k++) {
      maxi = NINFINITY; best = -1;     
      
      // On commence a coder (Start),si en phase ca vient d'une UTR 5' forward
      PICOMP((position % 3 == k),Start,Forward,UTR5F);
      // Il y a une insertion (frameshift). Saut de position de nucléotide ignoré.
      PICOMP(true,Ins,Forward,InitF1+(k+1)%3);
      // Il y a une deletion (frameshift)
      PICOMP(true,Del,Forward,InitF1+(k+2)%3);
    
      // On va tout droit.
      // S'il y  a un STOP en phase on ne peut continuer
      if ((position % 3 == k))
	LBP[InitF1+k].Update(Data.sig[DATA::Stop].weight[Signal::ForwardNo]);
      LBP[InitF1+k].Update(Data.sig[DATA::Don].weight[Signal::ForwardNo]);

      INSERT(InitF1+k);
    }
    // ----------------------------------------------------------------
    // ------------------------- Inits en reverse ---------------------
    // ----------------------------------------------------------------
    for (k = 3; k<6; k++) {
      maxi = NINFINITY; best = -1;

      // - on recommence a coder (Donneur) Ca vient d'un intron (no spliceable stop)
      PICOMPEN(true,Don,Reverse,IntronR1+((Data_Len-position-k) % 3),
	       ((((Data_Len-position) % 3) == k-3) ? Data.sig[DATA::Stop].weight[Signal::ReverseNo] : 0.0));

      // Not AfterG
      PICOMPEN(((Data_Len-position-k) % 3) == 2,Don,Reverse,IntronR3G ,
	       ((StartStop & DNASeq::isTAr) ? SplicedStopPen : 0.0));
      // Not AfterA
      PICOMPEN(((Data_Len-position-k) % 3) == 2,Don,Reverse,IntronR3A ,
	       ((StartStop & (DNASeq::isTGr | DNASeq::isTAr))  ? SplicedStopPen : 0.0));
      // Not AfterAG
      PICOMPEN(((Data_Len-position-k) % 3) == 1,Don,Reverse,IntronR2AG,
	       ((StartStop & DNASeq::isTr) ? SplicedStopPen : 0.0));

      // Il y a une insertion (frameshift)
      PICOMP(true,Ins,Reverse, InitR1+(k+2)%3);
      // Il y a une deletion (frameshift)
      PICOMP(true,Del,Reverse, InitR1+(k+1)%3);

      // On va tout droit.
      // S'il y  a un STOP en phase on ne peut continuer
      if (((Data_Len-position) % 3 == k-3)) 
	LBP[InitF1+k].Update(Data.sig[DATA::Stop].weight[Signal::ReverseNo]);
      LBP[InitF1+k].Update(Data.sig[DATA::Don].weight[Signal::ReverseNo]);

      INSERT(InitF1+k);
    }
    // ----------------------------------------------------------------
    // ------------------ Sngl en forward ----------------------------
    // ----------------------------------------------------------------
    for (k = 0; k < 3; k++) {
      maxi = NINFINITY; best = -1;     
      
      // On commence a coder (Start),si en phase ca vient d'une UTR 5' forward
      PICOMP((position % 3 == k),Start,Forward,UTR5F);
      // Il y a une insertion (frameshift). Saut de positionléotide ignore.
      PICOMP(true,Ins,Forward,SnglF1+(k+1)%3);
      // Il y a une deletion (frameshift)
      PICOMP(true,Del,Forward,SnglF1+(k+2)%3);
    
      // On va tout droit.
      // S'il y  a un STOP en phase on ne peut continuer
      if ((position % 3 == k))
	LBP[SnglF1+k].Update(Data.sig[DATA::Stop].weight[Signal::ForwardNo]);
      LBP[SnglF1+k].Update(Data.sig[DATA::Don].weight[Signal::ForwardNo]);
      
      INSERT(SnglF1+k);
    }
    // ----------------------------------------------------------------
    // ------------------------- Sngl en reverse ---------------------
    // ----------------------------------------------------------------
    for (k = 3; k<6; k++) {
      maxi = NINFINITY; best = -1;
      
      // On commence a coder (Stop). Ca vient d'une UTR 3' reverse
      PICOMP(((Data_Len-position) % 3 == k-3),Stop,Reverse,UTR3R);
      // Il y a une insertion (frameshift)
      PICOMP(true,Ins,Reverse, SnglR1+(k+2)%3);
      // Il y a une deletion (frameshift)
      PICOMP(true,Del,Reverse, SnglR1+(k+1)%3);

      // On va tout droit.
      // S'il y  a un STOP en phase on ne peut continuer
      if (((Data_Len-position) % 3 == k-3)) 
	LBP[SnglF1+k].Update(Data.sig[DATA::Stop].weight[Signal::ReverseNo]);
      LBP[SnglF1+k].Update(Data.sig[DATA::Don].weight[Signal::ReverseNo]);

      INSERT(SnglF1+k);
    }
    // ----------------------------------------------------------------
    // ------------------ Intrs en forward ----------------------------
    // ----------------------------------------------------------------
    for (k = 0; k < 3; k++) {
      maxi = NINFINITY; best = -1;     

      // On recommence a coder (Accepteur). Ca vient d'un intron (no spliceable stop)
      PICOMP(true,Acc,Forward,IntronF1+((position-k+3) % 3));

      // Not AfterT
      PICOMPEN(((position-k+3) % 3) == 1,Acc,Forward,IntronF2T ,
	       ((StopStop & (DNASeq::isGAf | DNASeq::isARf))  ? SplicedStopPen : 0.0));
      // Not AfterTG
      PICOMPEN(((position-k+3) % 3) == 2,Acc,Forward,IntronF3TG,
	       ((StopStop & DNASeq::isAf) ? SplicedStopPen : 0.0));
      // Not AfterTA
      PICOMPEN(((position-k+3) % 3) == 2,Acc,Forward,IntronF3TA,
	       ((StopStop & (DNASeq::isGf | DNASeq::isAf)) ? SplicedStopPen : 0.0));

      // Il y a une insertion (frameshift). Saut de position de nucléotide ignoré.
      PICOMP(true,Ins,Forward,IntrF1+(k+1)%3);
      // Il y a une deletion (frameshift)
      PICOMP(true,Del,Forward,IntrF1+(k+2)%3);
    
      // On va tout droit.
      // S'il y  a un STOP en phase on ne peut continuer
      if ((position % 3 == k))
	LBP[IntrF1+k].Update(Data.sig[DATA::Stop].weight[Signal::ForwardNo]);
      LBP[IntrF1+k].Update(Data.sig[DATA::Don].weight[Signal::ForwardNo]);

      INSERT(IntrF1+k);
    }
    // ----------------------------------------------------------------
    // ------------------------- Intrs en reverse ---------------------
    // ----------------------------------------------------------------
    for (k = 3; k<6; k++) {
      maxi = NINFINITY; best = -1;
      
      // - on recommence a coder (Donneur) Ca vient d'un intron (no spliceable stop)
      PICOMPEN(true,Don,Reverse,IntronR1+((Data_Len-position-k) % 3),((((Data_Len-position) % 3) == k-3) ? Data.sig[DATA::Stop].weight[Signal::ReverseNo] : 0.0));
      
      // Not AfterG
      PICOMPEN(((Data_Len-position-k) % 3) == 2,Don,Reverse,IntronR3G ,
	       ((StartStop & DNASeq::isTAr)  ? SplicedStopPen : 0.0));
      // Not AfterA
      PICOMPEN(((Data_Len-position-k) % 3) == 2,Don,Reverse,IntronR3A ,
	       ((StartStop & (DNASeq::isTAr | DNASeq::isTGr))  ? SplicedStopPen : 0.0));
      // Not AfterAG
      PICOMPEN(((Data_Len-position-k) % 3) == 1,Don,Reverse,IntronR2AG,
	       ((StartStop & DNASeq::isTr) ? SplicedStopPen : 0.0));

      // Il y a une insertion (frameshift)
      PICOMP(true,Ins,Reverse, IntrR1+(k+2)%3);
      // Il y a une deletion (frameshift)
      PICOMP(true,Del,Reverse, IntrR1+(k+1)%3);

      // On va tout droit.
      // S'il y  a un STOP en phase on ne peut continuer
      if (((Data_Len-position) % 3 == k-3)) 
	LBP[IntrF1+k].Update(Data.sig[DATA::Stop].weight[Signal::ReverseNo]);
      LBP[IntrF1+k].Update(Data.sig[DATA::Don].weight[Signal::ReverseNo]);

      INSERT(IntrF1+k);
    }
    // ----------------------------------------------------------------
    // ------------------ Terms en forward ----------------------------
    // ----------------------------------------------------------------
    for (k = 0; k < 3; k++) {
      maxi = NINFINITY; best = -1;     
      
      // On recommence a coder (Accepteur). Ca vient d'un intron
      PICOMP(true,Acc,Forward,IntronF1+((position-k+3) % 3));

      // Not AfterT
      PICOMPEN(((position-k+3) % 3) == 1,Acc,Forward,IntronF2T ,
	       ((StopStop & (DNASeq::isGAf | DNASeq::isARf))  ? SplicedStopPen :0.0));
      // Not AfterTG
      PICOMPEN(((position-k+3) % 3) == 2,Acc,Forward,IntronF3TG,
	       ((StopStop & DNASeq::isAf) ? SplicedStopPen : 0.0));
      // Not AfterTA
      PICOMPEN(((position-k+3) % 3) == 2,Acc,Forward,IntronF3TA,
	       ((StopStop & (DNASeq::isAf | DNASeq::isGf)) ? SplicedStopPen : 0.0));

      // Il y a une insertion (frameshift). Saut de positionléotide ignore.
      PICOMP(true,Ins,Forward,TermF1+(k+1)%3);
      // Il y a une deletion (frameshift)
      PICOMP(true,Del,Forward,TermF1+(k+2)%3);
    
      // On va tout droit.
      // S'il y  a un STOP en phase on ne peut continuer
      if ((position % 3 == k))
	LBP[TermF1+k].Update(Data.sig[DATA::Stop].weight[Signal::ForwardNo]);
      LBP[TermF1+k].Update(Data.sig[DATA::Don].weight[Signal::ForwardNo]);

      INSERT(TermF1+k);
    }
    // ----------------------------------------------------------------
    // ------------------------- Terms en reverse ---------------------
    // ----------------------------------------------------------------
    for (k = 3; k<6; k++) {
      maxi = NINFINITY; best = -1;
      
      // On commence a coder (Stop). Ca vient d'une UTR 3' reverse
      PICOMP(((Data_Len-position) % 3 == k-3),Stop,Reverse,UTR3R);
      // Il y a une insertion (frameshift)
      PICOMP(true,Ins,Reverse, TermR1+(k+2)%3);
      // Il y a une deletion (frameshift)
      PICOMP(true,Del,Reverse, TermR1+(k+1)%3);

      // On va tout droit.
      // S'il y  a un STOP en phase on ne peut continuer
      if (((Data_Len-position) % 3 == k-3)) 
	LBP[TermF1+k].Update(Data.sig[DATA::Stop].weight[Signal::ReverseNo]);
      LBP[TermF1+k].Update(Data.sig[DATA::Don].weight[Signal::ReverseNo]);

      INSERT(TermF1+k);
    }
    // ----------------------------------------------------------------
    // ------------------------ Intergenique --------------------------
    // ----------------------------------------------------------------
    // Ca peut venir d'une fin de 3' direct ou de 5' reverse
    maxi = NINFINITY; best = -1;
    
    // From 5' reverse
    PICOMP(true,tStart,Reverse, UTR5R);
    // From 3' direct
    PICOMP(true,tStop,Forward, UTR3F);
 
    // On reste intergenique
    // et les transstartNO/TransstopNO ???

    INSERT(InterGen);
    // ----------------------------------------------------------------
    // ---------------------- UTR 5' direct ---------------------------
    // ----------------------------------------------------------------
    maxi = NINFINITY; best = -1;
    
    // On vient de l'intergenique. 
    PICOMP(true,tStart,Forward, InterGen);
    // On peut venir aussi d'un intron d'UTR.
    PICOMP(true,Acc,Forward, IntronU5F);

    // On reste 5' direct. On ne prend pas le Start eventuel.
    LBP[UTR5F].Update(Data.sig[DATA::Start].weight[Signal::ForwardNo]);

    INSERT(UTR5F);
    // ----------------------------------------------------------------
    // ------------------- Introns d'UTR5F ----------------------------
    // ----------------------------------------------------------------
    maxi = NINFINITY; best = -1;
      
    // on quitte l'UTR 
    PICOMP(true,Don,Forward, UTR5F);

    // On reste intronique
    LBP[IntronU5F].Update(Data.sig[DATA::Acc].weight[Signal::ForwardNo]);

    INSERT(IntronU5F);
    // ----------------------------------------------------------------
    // ------------------- Introns d'UTR3F ----------------------------
    // ----------------------------------------------------------------
    maxi = NINFINITY; best = -1;
      
    // on quitte l'UTR 
    PICOMP(true,Don,Forward, UTR3F);

    // On reste intronique
    LBP[IntronU3F].Update(Data.sig[DATA::Acc].weight[Signal::ForwardNo]);

    INSERT(IntronU3F);
    // ----------------------------------------------------------------
    // ---------------------- UTR 3' direct ---------------------------
    // ----------------------------------------------------------------
    maxi = NINFINITY; best = -1;

    // Ca vient d'un Term ou d'un Sngl direct + STOP
    PICOMP(true,Stop,Forward, TermF1+position%3);
    PICOMP(true,Stop,Forward, SnglF1+position%3);
    // On peut venir aussi d'un intron d'UTR.
    PICOMP(true,Acc,Forward, IntronU3F);

    // On reste 3' direct

    INSERT(UTR3F);
    // ----------------------------------------------------------------
    // ----------------------- UTR 5'reverse --------------------------
    // ----------------------------------------------------------------
    maxi = NINFINITY; best = -1;
    
    // Ca vient d'un Sngl ou Init reverse + START
    PICOMP(true,Start,Reverse, InitR1+((Data_Len-position) % 3));
    PICOMP(true,Start,Reverse, SnglR1+((Data_Len-position) % 3));
    // On peut venir aussi d'un intron d'UTR.
    PICOMP(true,Don,Reverse, IntronU5R);

    // On reste 5' reverse
    LBP[UTR5R].Update(Data.sig[DATA::Start].weight[Signal::ReverseNo]);

    INSERT(UTR5R);
    // ----------------------------------------------------------------
    // ------------------- Introns d'UTR5R ----------------------------
    // ----------------------------------------------------------------
    maxi = NINFINITY; best = -1;
    
    // On quitte une UTR5R
    PICOMP(true,Acc,Reverse, UTR5R);

    // On reste intronique
    LBP[IntronU5R].Update(Data.sig[DATA::Acc].weight[Signal::ReverseNo]);

    INSERT(IntronU5R);
    // ----------------------------------------------------------------
    // ------------------- Introns d'UTR3R ----------------------------
    // ----------------------------------------------------------------
    maxi = NINFINITY; best = -1;
    
    // On quitte une UTR5R
    PICOMP(true,Acc,Reverse, UTR3R);

    // On reste intronique
    LBP[IntronU3R].Update(Data.sig[DATA::Acc].weight[Signal::ReverseNo]);

    INSERT(IntronU3R);
    // ----------------------------------------------------------------
    // ----------------------- UTR 3'reverse --------------------------
    // ----------------------------------------------------------------
    maxi = NINFINITY; best = -1;
    
    // On demarre depuis l'intergenique
    PICOMP(true,tStop,Reverse, InterGen);
    // On peut venir aussi d'un intron d'UTR.
    PICOMP(true,Don,Reverse, IntronU3R);

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
	PICOMPEN(true,Don,Forward, InitF1+((position-k+3) % 3),
		 (k == 0 ? Data.sig[DATA::Stop].weight[Signal::ForwardNo] : 0.0));
	PICOMPEN(true,Don,Forward, IntrF1+((position-k+3) % 3),
		 (k == 0 ? Data.sig[DATA::Stop].weight[Signal::ForwardNo] : 0.0));
      }

      // On reste intronique
      LBP[IntronF1+k].Update(Data.sig[DATA::Acc].weight[Signal::ForwardNo]);
      
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
      PICOMP(true,Don,Forward, InitF1+((position-k+3) % 3));
      PICOMP(true,Don,Forward, IntrF1+((position-k+3) % 3));
    }
    // On reste intronique
    LBP[IntronF2T].Update(Data.sig[DATA::Acc].weight[Signal::ForwardNo]);
    INSERT(IntronF2T);

    //
    // --- Intron Phase 2 after an TG|A
    //
    maxi = NINFINITY; best = -1;
    k = 2;
    if (StartStop & DNASeq::isTGf) {
      // - on quitte un Init ou un Intr
      PICOMP(true,Don,Forward, InitF1+((position-k+3) % 3));
      PICOMP(true,Don,Forward, IntrF1+((position-k+3) % 3));
    }
    // On reste intronique
    LBP[IntronF3TG].Update(Data.sig[DATA::Acc].weight[Signal::ForwardNo]);
    INSERT(IntronF3TG);

    //
    // --- Intron Phase 2 after a TA(A|G)
    //
    maxi = NINFINITY; best = -1;
    k = 2;
    if (StartStop & DNASeq::isTAf) {
      // - on quitte un Init ou un Intr
      PICOMP(true,Don,Forward, InitF1+((position-k+3) % 3));
      PICOMP(true,Don,Forward, IntrF1+((position-k+3) % 3));
    }
    // On reste intronique
    LBP[IntronF3TA].Update(Data.sig[DATA::Acc].weight[Signal::ForwardNo]);
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
	PICOMP(true,Acc,Reverse, IntrR1+((Data_Len-position-k) % 3));
	PICOMP(true,Acc,Reverse, TermR1+((Data_Len-position-k) % 3));
	}

      // On reste intronique
      LBP[IntronR1+k].Update(Data.sig[DATA::Acc].weight[Signal::ReverseNo]);

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
      PICOMP(true,Acc,Reverse, IntrR1+((Data_Len-position-k) % 3));
      PICOMP(true,Acc,Reverse, TermR1+((Data_Len-position-k) % 3));
    }
    // On reste intronique
    LBP[IntronR3G].Update(Data.sig[DATA::Acc].weight[Signal::ReverseNo]);
    INSERT(IntronR3G);
    
    //
    // --- Intron Phase 1 after an A (GT|AT)
    //
    maxi = NINFINITY; best = -1;
    k = 2;
    if (StopStop & DNASeq::isAr) {
      // - on quitte un Intr ou un Term
      PICOMP(true,Acc,Reverse, IntrR1+((Data_Len-position-k) % 3));
      PICOMP(true,Acc,Reverse, TermR1+((Data_Len-position-k) % 3));
    }
    // On reste intronique
    LBP[IntronR3A].Update(Data.sig[DATA::Acc].weight[Signal::ReverseNo]);
    INSERT(IntronR3A);

    //
    // --- Intron Phase 2 after an AG, AA ou GA (T)
    //
    maxi = NINFINITY; best = -1;
    k = 1;
if (StopStop & (DNASeq::isGAr | DNASeq::isARr)) {
      // - on quitte un Intr ou un Term
      PICOMP(true,Acc,Reverse, IntrR1+((Data_Len-position-k) % 3));
      PICOMP(true,Acc,Reverse, TermR1+((Data_Len-position-k) % 3));
    }
    // On reste intronique
    LBP[IntronR2AG].Update(Data.sig[DATA::Acc].weight[Signal::ReverseNo]);
    INSERT(IntronR2AG);

    // Account for Score and length penalty 
    ApplyScore(position, Data, 0);
    ApplyLengthPenalty(position, Data, 0);
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

  ApplyScore(position, Data, NoContentsUpdate);
  ApplyLengthPenalty(position, Data, NoContentsUpdate);
  
  // Avoid rounding errors on long sequences
  Normalize();
  
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
      INEED(IntronF1);
      INEED(IntronF2);
      INEED(IntronF3);
      INEED(IntronF2T);
      INEED(IntronF3TG);
      INEED(IntronF3TA);
    }
    // ---------------------------
    // 2- Exons F (from UTR3F)
    // ---------------------------
    if (ISPOSSIBLE(Stop,Forward))           
      INEED(UTR3F);

    // ---------------------------
    // 2b- Exons F (from frameshift)
    // ---------------------------
    if (ISPOSSIBLE(Ins,Forward) || ISPOSSIBLE(Del,Forward)) {
      INEED(InitF1); INEED(InitF2); INEED(InitF3);
      INEED(IntrF1); INEED(IntrF2); INEED(IntrF3);
      INEED(TermF1); INEED(TermF2); INEED(TermF3); 
      INEED(SnglF1); INEED(SnglF2); INEED(SnglF3); 
    }

    // ---------------------------
    // 3- Intron F (from ExonF)
    // ---------------------------
    if (ISPOSSIBLE(Acc,Forward)) {
      INEED(IntrF1);
      INEED(IntrF2);
      INEED(IntrF3);
      INEED(TermF1);
      INEED(TermF2);
      INEED(TermF3);
    }
    // ---------------------------
    // 4- UTR5F (from SnglF) - check
    // ---------------------------
    if (ISPOSSIBLE(Start,Forward)) {
      INEED(SnglF1+((position)%3));
      INEED(InitF1+((position)%3));
    }
    // ---------------------------
    // 5- UTR5F (from IntronU5F)
    // ---------------------------
    if (ISPOSSIBLE(Don,Forward))
      INEED(IntronU5F);
    // ---------------------------
    // 6- IntronU5F (from UTR5F)
    // ---------------------------
    if (ISPOSSIBLE(Acc,Forward))
      INEED(UTR5F);
    // ---------------------------
    // 7- UTR3F (from InterGen)
    // ---------------------------
    if (ISPOSSIBLE(tStop,Forward))
      INEED(InterGen);
    // ---------------------------
    // 8- UTR3F (from IntronU3F)
    // ---------------------------
    if (ISPOSSIBLE(Don,Forward))
      INEED(IntronU3F);
    // ---------------------------
    // 9- IntronU3F (from UTR3F)
    // ---------------------------
    if (ISPOSSIBLE(Acc,Forward))
      INEED(UTR3F);
    // ---------------------------
    // 10- InterG (from UTR5F)
    // ---------------------------
    if (ISPOSSIBLE(tStart,Forward))
      INEED(UTR5F);
    // ---------------------------
    // 11- InterG (from UTR3R)
    // ---------------------------
    if (ISPOSSIBLE(tStop,Reverse))
      INEED(UTR3R);
    // ---------------------------
    // 12- IntronU5R (from UTR5R)
    // ---------------------------
    if (ISPOSSIBLE(Don,Reverse))
      INEED(UTR5R);
    // ---------------------------
    // 13- UTR5R (from InterGen)
    // ---------------------------
    if (ISPOSSIBLE(tStart,Reverse))
      INEED(InterGen);
    // ---------------------------
    // 14- UTR5R (from IntronU5R)
    // ---------------------------
    if (ISPOSSIBLE(Acc,Reverse))
      INEED(IntronU5R);
    // ---------------------------
    // 15- IntronU3R (from UTR3R)
    // ---------------------------
    if (ISPOSSIBLE(Don,Reverse))
      INEED(UTR3R);
    // ---------------------------
    // 16- UTR3R (from IntronU3R)
    // ---------------------------
    if (ISPOSSIBLE(Acc,Reverse))
      INEED(IntronU3R);
    // ---------------------------
    // 17- UTR3R (from ExonR)
    // ---------------------------
    if (ISPOSSIBLE(Stop,Reverse)) {
      INEED(SnglR1+(Data_Len-position)%3);
      INEED(TermR1+(Data_Len-position)%3);
    }
    // ---------------------------
    // 18- ExonR (from UTR5R)
    // ---------------------------
    if (ISPOSSIBLE(Start,Reverse))
      INEED(UTR5R);
    // ---------------------------
    // 19- ExonR (from IntronR)
    // ---------------------------
    if (ISPOSSIBLE(Acc,Reverse)){
      INEED(IntronR1);
      INEED(IntronR2);
      INEED(IntronR3);
      INEED(IntronR2AG);
      INEED(IntronR3A);
      INEED(IntronR3G);
    }
    // ---------------------------
    // 19b- ExonR (from frameshift)
    // ---------------------------
    if (ISPOSSIBLE(Ins,Reverse) || ISPOSSIBLE(Del,Reverse)) {
      INEED(InitR1); INEED(InitR2); INEED(InitR3);
      INEED(IntrR1); INEED(IntrR2); INEED(IntrR3);
      INEED(TermR1); INEED(TermR2); INEED(TermR3); 
      INEED(SnglR1); INEED(SnglR2); INEED(SnglR3); 
    }

    // ---------------------------
    // 20- IntronR (from ExonR)
    // ---------------------------
    if (ISPOSSIBLE(Don,Reverse)) {
      INEED(InitR1);
      INEED(InitR2);
      INEED(InitR3);
      INEED(IntrR1);
      INEED(IntrR2);
      INEED(IntrR3);
    }

#ifdef DEBUGME
    for (k=0; k<(NbTracks-1)/2;k++)
      if (PrevBP[ForwardTracks[k]]) 
	printf("Cost %f Track %d a From %d\n",
	       PBest[ForwardTracks[k]],k,(PrevBP[ForwardTracks[k]]->State));

      if (PrevBP[UnorientedTracks[0]]) 
	printf("Cost %f Track %d From %d\n",
	       PBest[UnorientedTracks[0]],k,(PrevBP[UnorientedTracks[0]]->State));

    for (k=0; k<(NbTracks-1)/2;k++)
      if (PrevBP[ReverseTracks[k]]) 
	printf("Cost %f Track %d b From %d\n",
	       PBest[ReverseTracks[k]],k,(PrevBP[ReverseTracks[k]]->State));

    printf("---------- pos %d/%d, best %d norm %f ----------\n",Data_Len-position,position,(int)best,NormalizingPath);
#endif    

    // ----------------------------------------------------------------
    // ------------------------- Inits en reverse ---------------------
    // ----------------------------------------------------------------
    for (k = 0; k<3; k++) {
      maxi = NINFINITY; best = -1;
      
      // - on recommence a coder (Start) d'une UTR5R
      PICOMPR((Data_Len-position)%3 == k,Start,Reverse,UTR5R);
      // Il y a une insertion (frameshift). Saut de position de nucléotide ignoré.
      PICOMPR(true,Ins,Reverse,InitR1+(k+1)%3);
      // Il y a une deletion (frameshift)
      PICOMPR(true,Del,Reverse,InitR1+(k+2)%3);

      // On va tout droit.
      // S'il y  a un STOP en phase on ne peut continuer
      if (((Data_Len-position) % 3 == k)) 
	LBP[InitR1+k].Update(Data.sig[DATA::Stop].weight[Signal::ReverseNo]);
      LBP[InitR1+k].Update(Data.sig[DATA::Don].weight[Signal::ReverseNo]);

      INSERT(InitR1+k);
    }
    // ----------------------------------------------------------------
    // ------------------ Inits en forward ----------------------------
    // ----------------------------------------------------------------
    for (k = 0; k < 3; k++) {
      maxi = NINFINITY; best = -1;     
      
      // Début à partit d'un IntronF. Avoid Stop at end.
      PICOMPEN(true,Don,Forward,IntronF1+(position-k+3)%3,
	       ((position%3 == k) ? Data.sig[DATA::Stop].weight[Signal::ForwardNo] : 0.0));

      // Not after T (GA|AA|AG)
      PICOMPEN(((position-k+3)%3 == 1),Don,Forward,IntronF2T,
	       ((StartStop & DNASeq::isTf) ? SplicedStopPen : 0.0));
      // Not after TG|A
      PICOMPEN(((position-k+3)%3 == 2),Don,Forward,IntronF3TG,
	       ((StartStop & DNASeq::isTGf) ? SplicedStopPen : 0.0));
      // Not after TA(A|G)
      PICOMPEN(((position-k+3)%3 == 2),Don,Forward,IntronF3TA,
	       ((StartStop & DNASeq::isTAf) ? SplicedStopPen : 0.0));

      // Il y a une insertion (frameshift)
      PICOMPR(true,Ins,Forward, InitF1+(k+2)%3);
      // Il y a une deletion (frameshift)
      PICOMPR(true,Del,Forward, InitF1+(k+1)%3);

      // On va tout droit.
      // S'il y  a un STOP en phase on ne peut continuer
      if (((position)%3 == k))
	LBP[InitF1+k].Update(Data.sig[DATA::Stop].weight[Signal::ForwardNo]);
      LBP[InitF1+k].Update(Data.sig[DATA::Don].weight[Signal::ForwardNo]);

      INSERT(InitF1+k);
    }
    // ----------------------------------------------------------------
    // ------------------------- Sngl en reverse ---------------------
    // ----------------------------------------------------------------
    for (k = 0; k<3; k++) {
      maxi = NINFINITY; best = -1;
      
      // On commence a coder (Start). Ca vient d'une UTR 5' reverse
      PICOMPR(((Data_Len-position)%3 == k),Start,Reverse,UTR5R);
      // Il y a une insertion (frameshift). Saut de positionléotide ignore.
      PICOMPR(true,Ins,Reverse,SnglR1+(k+1)%3);
      // Il y a une deletion (frameshift)
      PICOMPR(true,Del,Reverse,SnglR1+(k+2)%3);
    
      // On va tout droit.
      // S'il y  a un STOP en phase on ne peut continuer
      if (((Data_Len-position)%3 == k)) 
	LBP[SnglR1+k].Update(Data.sig[DATA::Stop].weight[Signal::ReverseNo]);
      LBP[SnglR1+k].Update(Data.sig[DATA::Don].weight[Signal::ReverseNo]);

      INSERT(SnglR1+k);
    }
    // ----------------------------------------------------------------
    // ------------------ Sngl en forward ----------------------------
    // ----------------------------------------------------------------
    for (k = 0; k < 3; k++) {
      maxi = NINFINITY; best = -1;     
      
      // On commence a coder (Stop) sur une UTR3F
      PICOMPR(((position)%3 == k),Stop,Forward,UTR3F);
      // Il y a une insertion (frameshift)
      PICOMPR(true,Ins,Forward, SnglF1+(k+2)%3);
      // Il y a une deletion (frameshift)
      PICOMPR(true,Del,Forward, SnglF1+(k+1)%3);
      
      // On va tout droit.
      // S'il y  a un STOP en phase on ne peut continuer
      if (((position)%3 == k))
	LBP[SnglF1+k].Update(Data.sig[DATA::Stop].weight[Signal::ForwardNo]);
      LBP[SnglF1+k].Update(Data.sig[DATA::Don].weight[Signal::ForwardNo]);
      
      INSERT(SnglF1+k);
    }
    // ----------------------------------------------------------------
    // ------------------------- Intrs en reverse ---------------------
    // ----------------------------------------------------------------
    for (k = 0; k<3; k++) {
      maxi = NINFINITY; best = -1;
      
      // - on recommence a coder (Acc) Ca vient d'un intronR
      PICOMPR(true,Acc,Reverse,IntronR1+((Data_Len-position-k)%3));
      
      // Not AfterG
      PICOMPEN(((Data_Len-position-k)%3 == 2),Acc,Reverse,IntronR3G,
               ((StopStop & DNASeq::isGr)  ? SplicedStopPen : 0.0));
      // Not AfterA
      PICOMPEN(((Data_Len-position-k)%3 == 2),Acc,Reverse,IntronR3A,
               ((StopStop & DNASeq::isAr) ? SplicedStopPen : 0.0));
      // Not AfterAG
      PICOMPEN(((Data_Len-position-k)%3 == 1),Acc,Reverse,IntronR2AG,
               ((StopStop & (DNASeq::isGAr | DNASeq::isARr)) ? SplicedStopPen : 0.0));

      // Il y a une insertion (frameshift). Saut de position de nucléotide ignoré.
      PICOMPR(true,Ins,Reverse,IntrR1+(k+1)%3);
      // Il y a une deletion (frameshift)
      PICOMPR(true,Del,Reverse,IntrR1+(k+2)%3);
    
      // On va tout droit.
      // S'il y  a un STOP en phase on ne peut continuer
      if (((Data_Len-position)%3 == k)) 
	LBP[IntrR1+k].Update(Data.sig[DATA::Stop].weight[Signal::ReverseNo]);
      LBP[IntrR1+k].Update(Data.sig[DATA::Don].weight[Signal::ReverseNo]);

      INSERT(IntrR1+k);
    }
    // ----------------------------------------------------------------
    // ------------------ Intrs en forward ----------------------------
    // ----------------------------------------------------------------
    for (k = 0; k < 3; k++) {
      maxi = NINFINITY; best = -1;     

      // On recommence a coder (Donneur). Ca vient d'un intronF
      PICOMPEN(true,Don,Forward,IntronF1+((position-k+3) % 3),
	       ((position%3) == k) ? Data.sig[DATA::Stop].weight[Signal::ForwardNo] : 0.0);

      // Not AfterG
      PICOMPEN(((position-k+3)%3 == 1),Don,Forward,IntronF2T,
               ((StartStop & DNASeq::isTf)  ? SplicedStopPen : 0.0));
      // Not AfterA
      PICOMPEN(((position-k+3)%3 == 2),Don,Forward,IntronF3TG ,
               ((StartStop & DNASeq::isTGf) ? SplicedStopPen : 0.0));
      // Not AfterAG
      PICOMPEN(((position-k+3)%3 == 2),Don,Forward,IntronF3TA,
               ((StartStop & DNASeq::isTAf) ? SplicedStopPen : 0.0));

      // Il y a une insertion (frameshift)
      PICOMPR(true,Ins,Forward, IntrF1+(k+2)%3);
      // Il y a une deletion (frameshift)
      PICOMPR(true,Del,Forward, IntrF1+(k+1)%3);

      // On va tout droit.
      // S'il y  a un STOP en phase on ne peut continuer
      if (((position)%3 == k))
	LBP[IntrF1+k].Update(Data.sig[DATA::Stop].weight[Signal::ForwardNo]);
      LBP[IntrF1+k].Update(Data.sig[DATA::Don].weight[Signal::ForwardNo]);

      INSERT(IntrF1+k);
    }
    // ----------------------------------------------------------------
    // ------------------------- Terms en reverse ---------------------
    // ----------------------------------------------------------------
    for (k = 0; k<3; k++) {
      maxi = NINFINITY; best = -1;
      
      // On commence a coder (Acc). Ca vient d'un intronR
      PICOMPR(true,Acc,Reverse,IntronR1+(Data_Len-position-k)%3);

      // Not AfterG
      PICOMPEN(((Data_Len-position-k)%3 == 2),Acc,Reverse,IntronR3G,
               ((StopStop & DNASeq::isGr)  ? SplicedStopPen : 0.0));
      // Not AfterA
      PICOMPEN(((Data_Len-position-k)%3 == 2),Acc,Reverse,IntronR3A,
               ((StopStop & DNASeq::isAr) ? SplicedStopPen : 0.0));
      // Not AfterAG
      PICOMPEN(((Data_Len-position-k)%3 == 1),Acc,Reverse,IntronR2AG,
               ((StopStop & (DNASeq::isGAr | DNASeq::isARr)) ? SplicedStopPen : 0.0));

      // Il y a une insertion (frameshift). Saut de positionléotide ignore.
      PICOMPR(true,Ins,Reverse,TermR1+(k+1)%3);
      // Il y a une deletion (frameshift)
      PICOMPR(true,Del,Reverse,TermR1+(k+2)%3);
    
      // On va tout droit.
      // S'il y  a un STOP en phase on ne peut continuer
      if (((Data_Len-position)%3 == k)) 
	LBP[TermR1+k].Update(Data.sig[DATA::Stop].weight[Signal::ReverseNo]);
      LBP[TermR1+k].Update(Data.sig[DATA::Don].weight[Signal::ReverseNo]);

      INSERT(TermR1+k);
    }
    // ----------------------------------------------------------------
    // ------------------ Terms en forward ----------------------------
    // ----------------------------------------------------------------
    for (k = 0; k < 3; k++) {
      maxi = NINFINITY; best = -1;     
      
      // On recommence a coder (Stop). Ca vient d'une UTR3R
      PICOMPR((position%3 ==k),Stop,Forward,UTR3F);
      // Il y a une insertion (frameshift)
      PICOMPR(true,Ins,Forward, TermF1+(k+2)%3);
      // Il y a une deletion (frameshift)
      PICOMPR(true,Del,Forward, TermF1+(k+1)%3);

      // On va tout droit.
      // S'il y  a un STOP en phase on ne peut continuer
      if (((position)%3 == k))
	LBP[TermF1+k].Update(Data.sig[DATA::Stop].weight[Signal::ForwardNo]);
      LBP[TermF1+k].Update(Data.sig[DATA::Don].weight[Signal::ForwardNo]);

      INSERT(TermF1+k);
    }
    // ----------------------------------------------------------------
    // ------------------------ Intergenique --------------------------
    // ----------------------------------------------------------------
    // Ca peut venir d'une fin de 3' direct ou de 5' reverse
    maxi = NINFINITY; best = -1;
    
    // From 5' forward
    PICOMPR(true,tStart,Forward, UTR5F);
    // From 3' reverse
    PICOMPR(true,tStop,Reverse, UTR3R);
 
    // On reste intergenique
    // et les transstartNO/TransstopNO ???

    INSERT(InterGen);
    // ----------------------------------------------------------------
    // ----------------------- UTR 5'reverse --------------------------
    // ----------------------------------------------------------------
    maxi = NINFINITY; best = -1;
    
    // Ca vient d'interG (tStartR)
    PICOMPR(true,tStart,Reverse, InterGen);
    // On peut venir aussi d'un intron d'UTR.
    PICOMPR(true,Acc,Reverse, IntronU5R);

    // On reste 5' reverse
    LBP[UTR5R].Update(Data.sig[DATA::Start].weight[Signal::ReverseNo]);

    INSERT(UTR5R);
    // ----------------------------------------------------------------
    // ------------------- Introns d'UTR5R ----------------------------
    // ----------------------------------------------------------------
    maxi = NINFINITY; best = -1;
    
    // On quitte une UTR5R
    PICOMPR(true,Don,Reverse, UTR5R);

    // On reste intronique
    LBP[IntronU5R].Update(Data.sig[DATA::Acc].weight[Signal::ReverseNo]);

    INSERT(IntronU5R);
    // ----------------------------------------------------------------
    // ------------------- Introns d'UTR3R ----------------------------
    // ----------------------------------------------------------------
    maxi = NINFINITY; best = -1;
    
    // On quitte une UTR5R
    PICOMPR(true,Don,Reverse, UTR3R);

    // On reste intronique
    LBP[IntronU3R].Update(Data.sig[DATA::Acc].weight[Signal::ReverseNo]);

    INSERT(IntronU3R);
    // ----------------------------------------------------------------
    // ----------------------- UTR 3'reverse --------------------------
    // ----------------------------------------------------------------
    maxi = NINFINITY; best = -1;
    
    // On demarre depuis Sngl ou Term (Stop)
    PICOMPR(true,Stop,Reverse, TermR1+(Data_Len-position)%3);
    PICOMPR(true,Stop,Reverse, SnglR1+(Data_Len-position)%3);
    // On peut venir aussi d'un intron d'UTR.
    PICOMPR(true,Acc,Reverse, IntronU3R);

    // On reste 3' reverse

    INSERT(UTR3R);
    // ----------------------------------------------------------------
    // ---------------------- UTR 5' direct ---------------------------
    // ----------------------------------------------------------------
    maxi = NINFINITY; best = -1;
    
    // On vient d'un exon (SnglF ou InitF) sur un Start
    PICOMPR(true,Start,Forward, InitF1+(position)%3);
    PICOMPR(true,Start,Forward, SnglF1+(position)%3);
    // On peut venir aussi d'un intron d'UTR.
    PICOMPR(true,Don,Forward, IntronU5F);

    // On reste 5' direct. On ne prend pas le Start eventuel.
    LBP[UTR5F].Update(Data.sig[DATA::Start].weight[Signal::ForwardNo]);

    INSERT(UTR5F);
    // ----------------------------------------------------------------
    // ------------------- Introns d'UTR5F ----------------------------
    // ----------------------------------------------------------------
    maxi = NINFINITY; best = -1;
      
    // on quitte l'UTR 
    PICOMPR(true,Acc,Forward, UTR5F);

    // On reste intronique
    LBP[IntronU5F].Update(Data.sig[DATA::Acc].weight[Signal::ForwardNo]);

    INSERT(IntronU5F);
    // ----------------------------------------------------------------
    // ------------------- Introns d'UTR3F ----------------------------
    // ----------------------------------------------------------------
    maxi = NINFINITY; best = -1;
      
    // on quitte l'UTR 
    PICOMPR(true,Acc,Forward, UTR3F);

    // On reste intronique
    LBP[IntronU3F].Update(Data.sig[DATA::Acc].weight[Signal::ForwardNo]);

    INSERT(IntronU3F);
    // ----------------------------------------------------------------
    // ---------------------- UTR 3' direct ---------------------------
    // ----------------------------------------------------------------
    maxi = NINFINITY; best = -1;

    // Ca vient d'InterG (tStop)
    PICOMPR(true,tStop,Forward, InterGen);
    // On peut venir aussi d'un intron d'UTR.
    PICOMPR(true,Don,Forward, IntronU3F);

    // On reste 3' direct

    INSERT(UTR3F);
    // ----------------------------------------------------------------
    // ----------------- Introns de phase -k reverse ------------------
    // ----------------------------------------------------------------
    for (k = 0; k<3; k++) {
      maxi = NINFINITY; best = -1;
      
      // On quitte un Intr ou un Init reverse
      // no spliceable stop: 
      if (!(((StartStop & DNASeq::isTr) && k == 1) ||
            ((StartStop & (DNASeq::isTAr|DNASeq::isTGr)) && k == 2))) {
	PICOMPEN(true,Don,Reverse, InitR1+((Data_Len-position-k)%3),
		 (k == 0 ? Data.sig[DATA::Stop].weight[Signal::ReverseNo] : 0.0));
	PICOMPEN(true,Don,Reverse, IntrR1+((Data_Len-position-k)%3),
		 (k == 0 ? Data.sig[DATA::Stop].weight[Signal::ReverseNo] : 0.0));
      }

      // On reste intronique
      LBP[IntronR1+k].Update(Data.sig[DATA::Acc].weight[Signal::ReverseNo]);

      INSERT(IntronR1+k);
    }
    // ----------------------------------------------------------------
    // ----------------- Introns reverse speciaux ---------------------
    // ----------------------------------------------------------------
    //
    // --- Intron Phase 1 [G]---[AT]
    //
    maxi = NINFINITY; best = -1;
    k = 2;
    if (StartStop & DNASeq::isTAr) {
      // - on quitte un Init ou un Intr
      PICOMPR(true,Don,Reverse, InitR1+((Data_Len-position-k) % 3));
      PICOMPR(true,Don,Reverse, IntrR1+((Data_Len-position-k) % 3));
    }
    // On reste intronique
    LBP[IntronR3G].Update(Data.sig[DATA::Acc].weight[Signal::ReverseNo]);
    INSERT(IntronR3G);

    //
    // --- Intron Phase 1  [A]---[GT|AT]
    //
    maxi = NINFINITY; best = -1;
    k = 2;
    if (StartStop & (DNASeq::isTGr | DNASeq::isTAr)) {
      // - on quitte un Intr ou un Term
      PICOMPR(true,Don,Reverse, InitR1+((Data_Len-position-k) % 3));
      PICOMPR(true,Don,Reverse, IntrR1+((Data_Len-position-k) % 3));
    }
    // On reste intronique
    LBP[IntronR3A].Update(Data.sig[DATA::Acc].weight[Signal::ReverseNo]);
    INSERT(IntronR3A);

    //
    // --- Intron Phase 2 after an [AG|AA|GA]---[T]
    //
    maxi = NINFINITY; best = -1;
    k = 1;
    if (StartStop & DNASeq::isTr) {
      // - on quitte un Intr ou un Term
      PICOMPR(true,Don,Reverse, InitR1+((Data_Len-position-k) % 3));
      PICOMPR(true,Don,Reverse, IntrR1+((Data_Len-position-k) % 3));
    }
    // On reste intronique
    LBP[IntronR2AG].Update(Data.sig[DATA::Acc].weight[Signal::ReverseNo]);
    INSERT(IntronR2AG);
    // ----------------------------------------------------------------
    // ---------------- Introns de phase k forward --------------------
    // ----------------------------------------------------------------
    for (k = 0; k<3; k++) {
      maxi = NINFINITY; best = -1;
      
      // - on quitte un Term ou un Intr sur un Acc
      // no spliceable stop: 
      if (!(((StopStop & (DNASeq::isAf|DNASeq::isGf)) && k == 2) ||
            ((StopStop & (DNASeq::isGAf | DNASeq::isARf)) && k == 1)))  {
	PICOMPR(true,Acc,Forward,IntrF1+(position-k+3)%3);
	PICOMPR(true,Acc,Forward,TermF1+(position-k+3)%3);
      }

      // On reste intronique
      LBP[IntronF1+k].Update(Data.sig[DATA::Acc].weight[Signal::ForwardNo]);
      
      INSERT(IntronF1+k);
    }
    // ----------------------------------------------------------------
    // ----------------- Introns forward speciaux ---------------------
    // ----------------------------------------------------------------
    //
    // --- Intron Phase 2 [TG]---[A]
    //
    maxi = NINFINITY; best = -1;
    k = 2;
    if (StopStop & DNASeq::isAf) {
      // - on quitte un  Intr ou un Term
      PICOMPR(true,Acc,Forward, IntrF1+((position-k+3) % 3));
      PICOMPR(true,Acc,Forward, TermF1+((position-k+3) % 3));
    }
    // On reste intronique
    LBP[IntronF3TG].Update(Data.sig[DATA::Acc].weight[Signal::ForwardNo]);
    INSERT(IntronF3TG);

    //
    // --- Intron Phase 2 after a [TA]---[A|G]
    //
    maxi = NINFINITY; best = -1;
    k = 2;
    if (StopStop & (DNASeq::isAf | DNASeq::isGf)) {
      // - on quitte un Intr ou un Term
      PICOMPR(true,Acc,Forward, IntrF1+((position-k+3) % 3));
      PICOMPR(true,Acc,Forward, TermF1+((position-k+3) % 3));
    }
    // On reste intronique
    LBP[IntronF3TA].Update(Data.sig[DATA::Acc].weight[Signal::ForwardNo]);
    INSERT(IntronF3TA);

    //
    // --- Intron Phase 1  [T]---[GA|AA|AG]. 
    //
    maxi = NINFINITY; best = -1;
    k = 1;
    if (StopStop & (DNASeq::isGAf | DNASeq::isARf)) {
      // - on quitte un Intr ou un Term
      PICOMPR(true,Acc,Forward, IntrF1+((position-k+3) % 3));
      PICOMPR(true,Acc,Forward, TermF1+((position-k+3) % 3));
    }
    // On reste intronique
    LBP[IntronF2T].Update(Data.sig[DATA::Acc].weight[Signal::ForwardNo]);
    INSERT(IntronF2T);
}
