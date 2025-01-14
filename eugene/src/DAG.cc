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
#include <vector>


//#define DEBUGME 1
#include "Param.h"
extern Parameters PAR;

// extern vectors need to know the active tracks
extern vector<short int> activeTracks; 
extern vector<bool>      isActiveTrack;
extern std::map<int, std::vector<int> > bicodingStates;

//---------STATIC MEMBERS ---------------------------------
DNASeq *DAG::TheSeq = NULL;

double DAG::ExPrior 		= 0.0;
double DAG::InPrior 		= 0.0;
double DAG::IGPrior 		= 0.0;
double DAG::FivePrior 		= 0.0;
double DAG::ThreePrior 		= 0.0;
double DAG::IntronFivePrior 	= 0.0;
double DAG::RnaPrior            = 0.0;
double DAG::BiCodingPrior       = 0.0;
double DAG::UIRPrior            = 0.0;
double DAG::SplicedStopPen 	= 0.0;
int       DAG::estuse 		= 0;
double DAG::NormalizingPath	= 0.0;
double DAG::PBest[NbTracks];
BackPoint *DAG::PrevBP[NbTracks];
// ----------------------------------------------------------------
// Default constructor.
// ----------------------------------------------------------------
DAG :: DAG ()
{
  StartPosition   = 0;
  EndPosition     = 0;
  EvidenceName[0] = 0;
  pred            = NULL;

  for (int i=0; i < activeTracks.size(); i++)
    LBP[activeTracks[i]].Path.InitState(activeTracks[i], StartPosition);
}

// ----------------------------------------------------------------
// Constructor for AltPredicton
// Cond: Used only once for Main DAG
// NOT USED
// ----------------------------------------------------------------
DAG :: DAG (DNASeq* Seq, Parameters &PAR) 
{
    StartPosition   = 0;
    EndPosition     = 0;
    EvidenceName[0] = 0;
    pred            = NULL;

    for (int i=0; i < activeTracks.size(); i++) 
    {
        LBP[activeTracks[i]].Path.InitState(activeTracks[i], StartPosition);
    }
    DAG::TheSeq = Seq;
    DAG::SplicedStopPen = -PAR.getD("EuGene.SplicedStopPen");
    DAG::ExPrior        = PAR.getD("EuGene.ExonPrior");
    DAG::InPrior        = PAR.getD("EuGene.IntronPrior");
    DAG::IGPrior        = PAR.getD("EuGene.InterPrior"); 
    DAG::FivePrior      = PAR.getD("EuGene.FivePrimePrior");
    DAG::ThreePrior     = PAR.getD("EuGene.ThreePrimePrior");
    DAG::RnaPrior       = PAR.getD("EuGene.RnaPrior");
    DAG::BiCodingPrior  = PAR.getD("EuGene.BiCodingPrior");
    DAG::UIRPrior       = PAR.getD("EuGene.UIRPrior");
    DAG::estuse         = PAR.getI("Sensor.Est.use");
    DAG::IntronFivePrior = InPrior;
    DAG::NormalizingPath = 0.0;
}

// ----------------------------------------------------------------
// Default constructor. 
// Cond: Used only once for Main DAG !
// ----------------------------------------------------------------
DAG :: DAG (int start, int end, Parameters &PAR, DNASeq* Seq)
{
  StartPosition = start;
  EndPosition   = end;
  strcpy(EvidenceName,"Main");

  // Init all the active tracks
  for (int i = 0; i < activeTracks.size(); i++) { 
	LBP[activeTracks[i]].Path.InitState(activeTracks[i],StartPosition);
  }

  pred = NULL;
  DAG::TheSeq = Seq;

  DAG::SplicedStopPen = -PAR.getD("EuGene.SplicedStopPen");
  DAG::ExPrior        = PAR.getD("EuGene.ExonPrior");
  DAG::InPrior        = PAR.getD("EuGene.IntronPrior");
  DAG::IGPrior        = PAR.getD("EuGene.InterPrior"); 
  DAG::FivePrior      = PAR.getD("EuGene.FivePrimePrior");
  DAG::ThreePrior     = PAR.getD("EuGene.ThreePrimePrior");
  DAG::RnaPrior       = PAR.getD("EuGene.RnaPrior");
  DAG::BiCodingPrior  = PAR.getD("EuGene.BiCodingPrior");
  DAG::UIRPrior       = PAR.getD("EuGene.UIRPrior");
  DAG::IntronFivePrior = InPrior;
  DAG::estuse = PAR.getI("Sensor.Est.use");
  DAG::NormalizingPath = 0.0;
}

// ----------------------------------------------------------------
// Constructor from another DAG. 
// Plug the extremities of the new one on pos start and end on the other
// ----------------------------------------------------------------
DAG :: DAG (int start, int end, DAG *RefDag,char* name)
{
   StartPosition = start;
   EndPosition   = end;
   pred          = NULL;

	// Initialize the path for all the tracks
   for (int i = 0; i < NbTracks; i++) 
	LBP[i].Path.InitState(i, StartPosition);	
   strcpy(EvidenceName,name);

   // Create all the tracks
   for (int i = 0; i < NbTracks; i++) {
     LBP[i] = new Track(&RefDag->LBP[i]);
   }
}

// ----------------------------------------------------------------
//  Default destructor.
// ----------------------------------------------------------------
DAG :: ~DAG ()
{
  for  (int i = 0;  i < activeTracks.size();  i ++) LBP[activeTracks[i]].Zap();
}

// ----------------------------------------------------------------
//  Clean the activeTracks
// ----------------------------------------------------------------
void DAG :: Clean() 
{
    pred = NULL;
    
    for  (int i = 0;  i < activeTracks.size();  i ++) LBP[activeTracks[i]].Clean();
    DAG::NormalizingPath = 0.0;
}

// ----------------------------------------------------------------
//  Init the DAG at the start position
//  NOT USED
// ----------------------------------------------------------------
void DAG :: Init(int start, int end) 
{
    StartPosition = start;
    EndPosition   = end;
    // Init all the active tracks
    for (int i = 0; i < activeTracks.size(); i++) { 
        LBP[activeTracks[i]].Path.InitState(activeTracks[i],StartPosition);
    }
}

// ---------------------------------------------------------------------------
// Initial and terminal costs
// ---------------------------------------------------------------------------
void DAG :: WeightThePrior()
{
  // to marginally prefer Intr/Term over Init/Sngl on the borders
  // A better (operational) solution would be to use the good order
  // for finding optimum origin
  const double epsilon = 1e-6;

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
  
  // Special introns
  LBP[IntronF2T].Update(log(InPrior/13.0)/2.0);
  LBP[IntronF2A].Update(log(InPrior/13.0)/2.0);
  LBP[IntronF3TG].Update(log(InPrior/13.0)/2.0);
  LBP[IntronF3TA].Update(log(InPrior/13.0)/2.0);
  LBP[IntronF3TC].Update(log(InPrior/13.0)/2.0);
  LBP[IntronF3AG].Update(log(InPrior/13.0)/2.0);
  LBP[IntronR3G].Update(log(InPrior/13.0)/2.0);
  LBP[IntronR3A].Update(log(InPrior/13.0)/2.0);
  LBP[IntronR2GA].Update(log(InPrior/13.0)/2.0);
  LBP[IntronR2AG].Update(log(InPrior/13.0)/2.0);
  LBP[IntronR2AA].Update(log(InPrior/13.0)/2.0);
  LBP[IntronR2CA].Update(log(InPrior/13.0)/2.0);
  LBP[IntronR2GG].Update(log(InPrior/13.0)/2.0);

  // Sngl 
  LBP[SnglF1].Update(log(ExPrior/6.0)/2.0);
  LBP[SnglF2].Update(log(ExPrior/6.0)/2.0);
  LBP[SnglF3].Update(log(ExPrior/6.0)/2.0);
  LBP[SnglR1].Update(log(ExPrior/6.0)/2.0);
  LBP[SnglR2].Update(log(ExPrior/6.0)/2.0);
  LBP[SnglR3].Update(log(ExPrior/6.0)/2.0);

  // Intr 
  LBP[IntrF1].Update(log(ExPrior/6.0)/2.0+epsilon);
  LBP[IntrF2].Update(log(ExPrior/6.0)/2.0+epsilon);
  LBP[IntrF3].Update(log(ExPrior/6.0)/2.0+epsilon);
  LBP[IntrR1].Update(log(ExPrior/6.0)/2.0+epsilon);
  LBP[IntrR2].Update(log(ExPrior/6.0)/2.0+epsilon);
  LBP[IntrR3].Update(log(ExPrior/6.0)/2.0+epsilon);

  // Term
  LBP[TermF1].Update(log(ExPrior/6.0)/2.0+epsilon);
  LBP[TermF2].Update(log(ExPrior/6.0)/2.0+epsilon);
  LBP[TermF3].Update(log(ExPrior/6.0)/2.0+epsilon);
  LBP[TermR1].Update(log(ExPrior/6.0)/2.0+epsilon);
  LBP[TermR2].Update(log(ExPrior/6.0)/2.0+epsilon);
  LBP[TermR3].Update(log(ExPrior/6.0)/2.0+epsilon);

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

// Introns d'UTR5
  LBP[RnaF].Update(log(RnaPrior/2.0)/2.0);
  LBP[RnaR].Update(log(RnaPrior/2.0)/2.0);

  // BiCoding
  LBP[SnglF1F2].Update(log(BiCodingPrior/6.0)/2.0);
  LBP[SnglF1F3].Update(log(BiCodingPrior/6.0)/2.0);
  LBP[SnglF2F3].Update(log(BiCodingPrior/6.0)/2.0);
  LBP[SnglR1R2].Update(log(BiCodingPrior/6.0)/2.0);
  LBP[SnglR1R3].Update(log(BiCodingPrior/6.0)/2.0);
  LBP[SnglR2R3].Update(log(BiCodingPrior/6.0)/2.0);

  LBP[SnglF1R1].Update(log(BiCodingPrior/6.0)/2.0);
  LBP[SnglF1R2].Update(log(BiCodingPrior/6.0)/2.0);
  LBP[SnglF1R3].Update(log(BiCodingPrior/6.0)/2.0);
  LBP[SnglF2R1].Update(log(BiCodingPrior/6.0)/2.0);
  LBP[SnglF2R2].Update(log(BiCodingPrior/6.0)/2.0);
  LBP[SnglF2R3].Update(log(BiCodingPrior/6.0)/2.0);
  LBP[SnglF3R1].Update(log(BiCodingPrior/6.0)/2.0);
  LBP[SnglF3R2].Update(log(BiCodingPrior/6.0)/2.0);
  LBP[SnglF3R3].Update(log(BiCodingPrior/6.0)/2.0);

  
  // UIR
  LBP[UIRF].Update(log(UIRPrior/2.0)/2.0);
  LBP[UIRR].Update(log(UIRPrior/2.0)/2.0); 
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
  LBP[IntronF2A].LoadPenalty(PAR.getC("EuGene.IntronDist"));
  LBP[IntronF3TG].LoadPenalty(PAR.getC("EuGene.IntronDist"));
  LBP[IntronF3TA].LoadPenalty(PAR.getC("EuGene.IntronDist"));
  LBP[IntronF3TC].LoadPenalty(PAR.getC("EuGene.IntronDist"));
  LBP[IntronF3AG].LoadPenalty(PAR.getC("EuGene.IntronDist"));
  LBP[IntronR3G].LoadPenalty(PAR.getC("EuGene.IntronDist"));
  LBP[IntronR3A].LoadPenalty(PAR.getC("EuGene.IntronDist"));
  LBP[IntronR2GA].LoadPenalty(PAR.getC("EuGene.IntronDist"));
  LBP[IntronR2AG].LoadPenalty(PAR.getC("EuGene.IntronDist"));
  LBP[IntronR2AA].LoadPenalty(PAR.getC("EuGene.IntronDist"));
  LBP[IntronR2CA].LoadPenalty(PAR.getC("EuGene.IntronDist"));
  LBP[IntronR2GG].LoadPenalty(PAR.getC("EuGene.IntronDist"));

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

  // Rna 
  LBP[RnaF].LoadPenalty(PAR.getC("EuGene.RnaDist"));
  LBP[RnaR].LoadPenalty(PAR.getC("EuGene.RnaDist"));

  // Bicoding 
  LBP[SnglF1F2].LoadPenalty(PAR.getC("EuGene.OverlapDist"));
  LBP[SnglF1F3].LoadPenalty(PAR.getC("EuGene.OverlapDist"));
  LBP[SnglF2F3].LoadPenalty(PAR.getC("EuGene.OverlapDist"));
  LBP[SnglR1R2].LoadPenalty(PAR.getC("EuGene.OverlapDist"));
  LBP[SnglR1R3].LoadPenalty(PAR.getC("EuGene.OverlapDist"));
  LBP[SnglR2R3].LoadPenalty(PAR.getC("EuGene.OverlapDist"));

  LBP[SnglF1R1].LoadPenalty(PAR.getC("EuGene.OverlapDist"));
  LBP[SnglF1R2].LoadPenalty(PAR.getC("EuGene.OverlapDist"));
  LBP[SnglF1R3].LoadPenalty(PAR.getC("EuGene.OverlapDist"));
  LBP[SnglF2R1].LoadPenalty(PAR.getC("EuGene.OverlapDist"));
  LBP[SnglF2R2].LoadPenalty(PAR.getC("EuGene.OverlapDist"));
  LBP[SnglF2R3].LoadPenalty(PAR.getC("EuGene.OverlapDist"));
  LBP[SnglF3R1].LoadPenalty(PAR.getC("EuGene.OverlapDist"));
  LBP[SnglF3R2].LoadPenalty(PAR.getC("EuGene.OverlapDist"));
  LBP[SnglF3R3].LoadPenalty(PAR.getC("EuGene.OverlapDist"));


  // UIR
  LBP[UIRF].LoadPenalty(PAR.getC("EuGene.UIRDist"));
  LBP[UIRR].LoadPenalty(PAR.getC("EuGene.UIRDist"));
  

}
// ----------------------------------------------------------------
//  Build the prediction by backtracing.
// ----------------------------------------------------------------
double DAG :: BuildPrediction (int From, int To, int Forward)
{
  double maxi,BestU,PBest[NbTracks];
  BackPoint *PrevBP[NbTracks];

  int j, k = (Forward ? To+2 :From-1); // SeqLen+1 and -1 typically
  vector<short int>::iterator it;

  // Insert best possible backpoint at the start of the algo (the
  // insert is not automatically possible, cf cost dist. 
  // on length)
  for ( it = activeTracks.begin(); it != activeTracks.end(); ++it) 
  {
    PrevBP[*it] = LBP[*it].BestUsable(k, &PBest[*it], 0);
    LBP[*it].ForceNew(*it, k, PBest[*it], PrevBP[*it]);
  }

  // Select where Backtrace should start
  j    = activeTracks[0];        // will contain the number of the tracks where start the backtrace
  maxi = PBest[activeTracks[0]]; // will contain the corresponding score
  for ( it = activeTracks.begin(); it != activeTracks.end(); ++it) 
  { 
    BestU = PBest[*it];
    // Un test tordu pour casser le cou aux NaN
    if (isnan(maxi) || (BestU > maxi)) {
      maxi = BestU;
      j    = *it;
    }
  }

  pred = LBP[j].BackTrace(From,To,Forward);
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
void DAG :: StatActive() 
{
  int BPalloc = 0;

  for (int k = 0; k < activeTracks.size(); k++) {
	BPalloc += LBP[activeTracks[k]].NumBPAlloc;
  }
  
   printf("(%d active BP, ",BPalloc);
}
// ----------------------------------------------------------------
//  Print collected BP stats
// ----------------------------------------------------------------
void DAG :: StatGC() 
{
  int BPcollect = 0;

  for (int k = 0; k < activeTracks.size(); k++) {
    BPcollect += LBP[activeTracks[k]].NumBPCollect;
  }
   printf("%d collected).\n",BPcollect);
}
// ----------------------------------------------------------------
//  MarkandSweep garbage collector
// ----------------------------------------------------------------
void DAG :: MarkAndSweep(int pos,int verbose, int latency)
{
  
  int k;
  int Horizon = pos;

  if (verbose) {
    printf("GC started ");
    StatActive();
  }

  // We first compute the maximum horizon that will be GC'd. This is
  // a purely approximative approach. There is a priori no guarantee
  // that everything will be GC'd with this but most things should.
  // For each track, the penalty distribution may have a different
  // length and everything inside this length is a source pointeur for
  // tracing. The horizon is therefore the Minimum over all Tracks, of
  // the maximum of 0 and (currentpos - PenD.MaxLen - k*GCLATENCY). k
  // = 2 should be nice.

  for (k = 0; k < activeTracks.size(); k++) 
    Horizon = Min(Horizon,Max(0,pos - LBP[activeTracks[k]].PenD->MaxLen - 2*latency));

  for (k = 0; k < activeTracks.size() ; k++)
    LBP[activeTracks[k]].ClearMark(Horizon);

  for (k = 0; k < activeTracks.size() ; k++) 
    LBP[activeTracks[k]].Mark(Horizon);
  
  for (k = 0; k < activeTracks.size() ; k++) 
    LBP[activeTracks[k]].Sweep(Horizon);
  

  if (verbose) 
    StatGC();
}

// ----------------------------------------------------------------
//  Normalize to avoid rounding errors
// ----------------------------------------------------------------
inline void DAG::Normalize()
{
  double BestU,maxi = -NINFINITY;

  for (int k = 0 ; k < activeTracks.size(); k++) {
    //fprintf(stderr, "Optimal : %f\n", LBP[activeTracks[k]].Optimal);
    BestU = LBP[activeTracks[k]].Optimal;
    if ((BestU > NINFINITY) && (BestU < maxi)) 
      maxi = BestU;
  }
  
  for (int k = 0 ; k < activeTracks.size(); k++) 
    LBP[activeTracks[k]].Update(-maxi);
  NormalizingPath += maxi;
  
}

// ----------------------------------------------------------------
//  Account for length penalty distributions
// ----------------------------------------------------------------
inline void DAG::ApplyLengthPenalty(int position, DATA Data, int NoContentsUpdate)
{
  int Data_Len = TheSeq->SeqLen;

  if ((position < Data_Len) && (!NoContentsUpdate))
    for (int k=0; k < activeTracks.size(); k++)
      LBP[activeTracks[k]].PayTheSlope();
}

// ----------------------------------------------------------------
//  Apply the supermix compute of the two scores
// ----------------------------------------------------------------
inline double DAG::SuperMix(double a, double b)
{
  double SuperP = 0.5;
  return SuperP*Max(a,b)+(1.0-SuperP)*Min(a,b);
}

// ----------------------------------------------------------------
//  Account for contents score
// ----------------------------------------------------------------
void DAG::ApplyScore ( int position, DATA Data, int NoContentsUpdate )
{
	int Data_Len = TheSeq->SeqLen;
	vector<short int>::iterator it;

	if ( ( position < Data_Len ) && ( !NoContentsUpdate ) )
	{
		for ( it = activeTracks.begin(); it != activeTracks.end(); ++it) 
		{
			if ( !State (*it).IsBicoding() )
			{
				LBP[*it].Update ( Data.contents[SensorContents[*it]] );
			}
			else // case where the track is bicoding
			{
				LBP[*it].Update (SuperMix (
					Data.contents[SensorContents[bicodingStates[*it][0]]],
					Data.contents[SensorContents[bicodingStates[*it][1]]]) );
			}
		}
	}
}
// ----------------------------------------------------------------
// A first set of macros
// ----------------------------------------------------------------
#define ISPOSSIBLE(X,Y) (!(isinf(Data.sig[DATA::X].weight[Y])))

#define INEED(K) if (!PrevBP[Strand ?  ReverseIt[K] : K] && isActiveTrack[K])		\
    (PrevBP[Strand ?  ReverseIt[K] : K] =				\
     LBP[Strand ? ReverseIt[K] : K].BestUsable(position, &PBest[Strand ? ReverseIt[K] : K]))

// ----------------------------------------------------------------
// Calcul des meilleures opening edges
// ----------------------------------------------------------------
inline void DAG::ComputeRequired(enum Signal::Edge Strand, DATA Data, int position)
{
  int Data_Len = TheSeq->SeqLen;
  int FrameF = (Strand ? ((Data_Len-position+3) % 3) : (position % 3));
  int FrameR = (Strand ? (position % 3) : ((Data_Len-position+3) % 3));

  enum Signal::Edge FStrand = (Strand ? Signal::Reverse : Signal::Forward);
  enum Signal::Edge RStrand = (Strand ? Signal::Forward : Signal::Reverse);

  for (int k = 0; k < NbTracks; k++) {
    PrevBP[k] = NULL;
    PBest[k]  = NINFINITY;
  }
  
  // ---------------------------
  // ExonsR spliced (from Intron)
  // ---------------------------
  if (ISPOSSIBLE(Don,RStrand)) {
    INEED(IntronR1); 
    INEED(IntronR2); 
    INEED(IntronR3);
    INEED(IntronR3G);
    INEED(IntronR3A);
    INEED(IntronR2GA);
    INEED(IntronR2AG);
    INEED(IntronR2AA);
    INEED(IntronR2CA);
    INEED(IntronR2GG); 
  }
  
  // ---------------------------
  // ExonsR (from UTR3R)
  // ---------------------------
  if (ISPOSSIBLE(Stop,RStrand)) 
    INEED(UTR3R);
    
  // ---------------------------
  // ExonsR (from frameshift)
  // ---------------------------
    if (ISPOSSIBLE(Ins1,RStrand) || ISPOSSIBLE(Ins2,RStrand) || ISPOSSIBLE(Ins3,RStrand) || 
        ISPOSSIBLE(Del1,RStrand) || ISPOSSIBLE(Del2,RStrand) || ISPOSSIBLE(Del3,RStrand) ) {
    INEED(InitR1); INEED(InitR2); INEED(InitR3);
    INEED(IntrR1); INEED(IntrR2); INEED(IntrR3);
    INEED(TermR1); INEED(TermR2); INEED(TermR3); 
    INEED(SnglR1); INEED(SnglR2); INEED(SnglR3); 
  }
  
  // ---------------------------
  // IntronR (from ExonR)
  // ---------------------------
  if (ISPOSSIBLE(Acc,RStrand)) {
    INEED(IntrR1); INEED(IntrR2); INEED(IntrR3);
    INEED(TermR1); INEED(TermR2); INEED(TermR3); 
  }

  // ---------------------------
  // UTR5R (from Init/SnglR) - check or UIRR (from SnglR)
  // ---------------------------
  if (ISPOSSIBLE(Start,RStrand)) {
    INEED(InitR1+FrameR);
    INEED(SnglR1+FrameR); 
  }

  // ---------------------------
  // UTR5R (from IntronU5R)
  // ---------------------------
  if (ISPOSSIBLE(Don,RStrand))
    INEED(IntronU5R);

  // ---------------------------
  // IntronU5R (from UTR5R)
  // ---------------------------
  if (ISPOSSIBLE(Acc,RStrand))
    INEED(UTR5R);

  // ---------------------------
  // UTR3R (from InterGen)
  // ---------------------------
  if (ISPOSSIBLE(tStop,RStrand))
    INEED(InterGen);

  // ---------------------------
  // UTR3R (from IntronU3R)
  // ---------------------------
  if (ISPOSSIBLE(Don,RStrand))
    INEED(IntronU3R);

  // ---------------------------
  // IntronU3R (from UTR3R)
  // ---------------------------
  if (ISPOSSIBLE(Acc,RStrand))
    INEED(UTR3R);

  // ---------------------------
  // InterG (from UTR5R)
  // ---------------------------
  if (ISPOSSIBLE(tStart,RStrand))
  {
    INEED(UTR5R);
  }

  // ---------------------------
  // InterG (from UTR3F)
  // ---------------------------
  if (ISPOSSIBLE(tStop,FStrand))
  {
    INEED(UTR3F);
  }

  // ---------------------------
  // InterG or UIRR (from RnaR)
  // ---------------------------
  if (ISPOSSIBLE(tStartNpc,RStrand))
  {
    INEED(RnaR);
  }
  // ---------------------------
  //  InterG or UIRF (from RnaF)
  // ---------------------------
  if (ISPOSSIBLE(tStopNpc,FStrand))
  {
    INEED(RnaF);
  }

  // ---------------------------
  // RnaF (from InterGen or UIRF)
  // ---------------------------
  if (ISPOSSIBLE(tStartNpc, FStrand))
  {
        INEED(InterGen);
        INEED(UIRF);
  }

   // ---------------------------
  // RnaR (from InterGen or UIRR)
  // ---------------------------
   if (ISPOSSIBLE(tStopNpc, RStrand))
   {
        INEED(InterGen);
        INEED(UIRR);
   }

  // ---------------------------
  // IntronU5F (from UTR5F)
  // ---------------------------
  if (ISPOSSIBLE(Don,FStrand))
    INEED(UTR5F);

  // ---------------------------
  // UTR5F (from InterGen)
  // ---------------------------
  if (ISPOSSIBLE(tStart,FStrand))
    INEED(InterGen);

  // ---------------------------
  // UTR5F (from IntronU5F)
  // ---------------------------
  if (ISPOSSIBLE(Acc,FStrand))
    INEED(IntronU5F);

  // ---------------------------
  // IntronU3R (from UTR3F)
  // ---------------------------
  if (ISPOSSIBLE(Don,FStrand))
    INEED(UTR3F);

  // ---------------------------
  // UTR3F (from IntronU3F)
  // ---------------------------
  if (ISPOSSIBLE(Acc,FStrand))
    INEED(IntronU3F);

  // ---------------------------
  // UTR3F (from ExonF) or UIRF (from SnglF)
  // ---------------------------
  if (ISPOSSIBLE(Stop,FStrand)) {
    INEED(TermF1+FrameF);
    INEED(SnglF1+FrameF); 
  }

  // ---------------------------
  // ExonF (from UTR5F)
  // ---------------------------
  if (ISPOSSIBLE(Start,FStrand))
    INEED(UTR5F);

  // ---------------------------
  // ExonF (from IntronF)
  // ---------------------------
  if (ISPOSSIBLE(Acc,FStrand)) {
    INEED(IntronF2T);
    INEED(IntronF2A);
    INEED(IntronF3TG);
    INEED(IntronF3TA);
    INEED(IntronF3TC);
    INEED(IntronF3AG); 
    INEED(IntronF1); 
    INEED(IntronF2); 
    INEED(IntronF3);
  }

  // ---------------------------
  // ExonF (from frameshift)
  // ---------------------------
   if (ISPOSSIBLE(Ins1,FStrand) || ISPOSSIBLE(Ins2,FStrand) || ISPOSSIBLE(Ins3,FStrand) || 
       ISPOSSIBLE(Del1,FStrand) || ISPOSSIBLE(Del2,FStrand) || ISPOSSIBLE(Del3,FStrand)) {
    INEED(InitF1); INEED(InitF2); INEED(InitF3);
    INEED(IntrF1); INEED(IntrF2); INEED(IntrF3);
    INEED(TermF1); INEED(TermF2); INEED(TermF3); 
    INEED(SnglF1); INEED(SnglF2); INEED(SnglF3); 
  }
  
  // ---------------------------
  // SnglRR (from SnglR)
  //----------------------------
  if (ISPOSSIBLE(Stop,RStrand)) {
	if (FrameR != 0) 
	INEED(SnglR1); 
        if (FrameR != 1)
	INEED(SnglR2);
        if (FrameR != 2)
	INEED(SnglR3);
 }
  // ---------------------------
  // SnglFF (from SnglF)
  // ---------------------------
  if (ISPOSSIBLE(Start,FStrand)) {
    if (FrameF != 0) 
	INEED(SnglF1); 
    if (FrameF != 1)
	INEED(SnglF2);
    if (FrameF != 2)
	INEED(SnglF3);
  }
  // ---------------------------
  // SnglFR (from SnglF)
  // ---------------------------
  if (ISPOSSIBLE(Stop, RStrand)) {
    INEED(SnglF1); INEED(SnglF2);INEED(SnglF3);	
  }

  // ---------------------------
  // SnglFR (from SnglR)
  // ---------------------------
  if (ISPOSSIBLE(Start, FStrand)) {
    INEED(SnglR1); INEED(SnglR2);INEED(SnglR3);	
  }

  // ---------------------------
  // IntronF (from ExonF)
  // ---------------------------
  if (ISPOSSIBLE(Don,FStrand)) {
    INEED(InitF1); INEED(InitF2); INEED(InitF3);
    INEED(IntrF1); INEED(IntrF2); INEED(IntrF3);
  }
  
 // ---------------------------
  // SnglF (from SnglFF)
  //----------------------------
  if (ISPOSSIBLE(Stop,FStrand)) {
	if (FrameF == 0 || (FrameF == 1)) 
		INEED(SnglF1F2);
	if (FrameF == 0 || (FrameF == 2))
		INEED(SnglF1F3);
    if (FrameF == 1 || (FrameF == 2))
		INEED(SnglF2F3);
 }

  // ---------------------------
  // SnglR (from SnglRR)
  // ---------------------------
  if (ISPOSSIBLE(Start,RStrand)) {
    if (FrameR == 0 || (FrameR == 1)) 
		INEED(SnglR1R2);
	if (FrameR == 0 || (FrameR == 2))
		INEED(SnglR1R3);
    if (FrameR == 1 || (FrameR == 2))
		INEED(SnglR2R3);
  }

  // ---------------------------
  // SnglF (from SnglFR)
  //----------------------------
  if (ISPOSSIBLE(Start,RStrand)) {
	if (FrameR == 0)
	{
		INEED(SnglF1R1);
		INEED(SnglF2R1);
		INEED(SnglF3R1);
	}
	else if (FrameR == 1)
	{ 
		INEED(SnglF1R2);
		INEED(SnglF2R2);
		INEED(SnglF3R2);
	}
	else if (FrameR == 2)
	{
		INEED(SnglF1R3);
		INEED(SnglF2R3);
		INEED(SnglF3R3);
	}
  }

  // ---------------------------
  // SnglR (from SnglFR)
  //----------------------------
  if (ISPOSSIBLE(Stop, FStrand)) {
    if (FrameF == 0)
	{
		INEED(SnglF1R1);
		INEED(SnglF1R2);
		INEED(SnglF1R3);
 	}
	else if (FrameF == 1)
	{
		INEED(SnglF2R1);
		INEED(SnglF2R2);
		INEED(SnglF2R3);
	}
	else if (FrameF == 2)
	{
		INEED(SnglF3R1);
		INEED(SnglF3R2);
		INEED(SnglF3R3);
	}
  }	

  // ---------------------------
  // SnglF (from UIRF)
  // ---------------------------
  if (ISPOSSIBLE(Start, FStrand)) {
    INEED(UIRF);
  }
  // ---------------------------
  // SnglR (from UIRR)
  // ---------------------------
  if (ISPOSSIBLE(Stop, RStrand)) {
    INEED(UIRR);
  }


#ifdef DEBUGME
    for (int k=0; k<(NbTracks-1)/2;k++)
      if (PrevBP[ForwardTracks[k]]) 
	printf("Cost %f Track %d a From %d\n",
	       PBest[ForwardTracks[k]],k,(PrevBP[ForwardTracks[k]]->State));

      if (PrevBP[UnorientedTracks[0]]) 
	printf("Cost %f Track %d From %d\n",
	       PBest[UnorientedTracks[0]],0,(PrevBP[UnorientedTracks[0]]->State));

    for (int k=0; k<(NbTracks-1)/2;k++)
      if (PrevBP[ReverseTracks[k]]) 
	printf("Cost %f Track %d b From %d\n",
	       PBest[ReverseTracks[k]],k,(PrevBP[ReverseTracks[k]]->State));

    printf("---------- pos %d/%d, norm %f ----------\n",position,Data_Len-position,NormalizingPath);
#endif
}

// ----------------------------------------------------------------
// Spliceable stops... there is probably a nice small automata for this
// ----------------------------------------------------------------
#ifndef TRANSL_TABLE_6
/*inline bool spliceStopPreTAx(int StartStop, enum Signal::Edge Strand)  //TAx
{
  if (Strand == Signal::Forward) 
    return (StartStop & DNASeq::isTAf);
  else
    return (StartStop & DNASeq::isTAr);
}
// ----------------------------------------------------------------
inline bool spliceStopPostTAx(int StopStop, enum Signal::Edge Strand)  //xx(G|A) xxG
{
  if (Strand == Signal::Forward) 
    return (StopStop & (DNASeq::isGf | DNASeq::isAf));
  else
    return (StopStop & DNASeq::isGr);
}
// ----------------------------------------------------------------
inline bool spliceStopPrexxA(int StartStop, enum Signal::Edge Strand)  //TGx (TG|TA)x
{
  if (Strand == Signal::Forward) 
    return (StartStop & DNASeq::isTGf);
  else
    return (StartStop & (DNASeq::isTGr|DNASeq::isTAr));
}
// ----------------------------------------------------------------
inline bool spliceStopPostxxA(int StopStop, enum Signal::Edge Strand)  //xxA
{
  if (Strand == Signal::Forward) 
    return (StopStop & DNASeq::isAf);
  else
    return (StopStop & DNASeq::isAr);
}
// ----------------------------------------------------------------
inline bool spliceStopPreTxx(int StartStop, enum Signal::Edge Strand)  //Txx
{
  if (Strand == Signal::Forward) 
    return (StartStop & DNASeq::isTf);
  else
    return (StartStop & DNASeq::isTr);
}
// ----------------------------------------------------------------
inline bool spliceStopPostTxx(int StopStop, enum Signal::Edge Strand)  //x(GA|A(A|G))
{
  if (Strand == Signal::Forward) 
    return (StopStop & (DNASeq::isGAf | DNASeq::isARf));
  else
    return (StopStop & (DNASeq::isGAr | DNASeq::isARr));
}*/
#endif

#ifdef TRANSL_TABLE_6
/*
// ----------------------------------------------------------------
// Only one stop codon TGA
// ----------------------------------------------------------------
inline bool spliceStopPreTAx(int StartStop, enum Signal::Edge Strand)  //TAx --> not a stop codon
{
	return false;
}
// ----------------------------------------------------------------
inline bool spliceStopPostTAx(int StopStop, enum Signal::Edge Strand)  //xx(G|A) xxG --> not a stop codon
{
  return false;
}
// ----------------------------------------------------------------
inline bool spliceStopPrexxA(int StartStop, enum Signal::Edge Strand)  //TGx (TG)x
{
  if (Strand == Signal::Forward) 
    return (StartStop & DNASeq::isTGf);
  else
    return (StartStop & DNASeq::isTGr);
}
// ----------------------------------------------------------------
inline bool spliceStopPostxxA(int StopStop, enum Signal::Edge Strand)  //xxA
{
  if (Strand == Signal::Forward) 
    return (StopStop & DNASeq::isAf);
  else
    return (StopStop & DNASeq::isAr);
}
// ----------------------------------------------------------------
inline bool spliceStopPreTxx(int StartStop, enum Signal::Edge Strand)  //Txx
{
  if (Strand == Signal::Forward) 
    return (StartStop & DNASeq::isTf);
  else
    return (StartStop & DNASeq::isTr);
}
// ----------------------------------------------------------------
inline bool spliceStopPostTxx(int StopStop, enum Signal::Edge Strand)  //x(GA)
{
  if (Strand == Signal::Forward) 
    return (StopStop & DNASeq::isGAf);
  else
    return (StopStop & DNASeq::isGAr);
}*/
#endif



inline bool beginSpliceStopF2(int StartStop)
{
  if  ( (StartStop & DNASeq::startTf) || (StartStop & DNASeq::startAf) )
     return true;
  return false;
}

inline bool beginSpliceStopF2(int StartStop, int state)
{
  if  ( ( state == IntronF2T && (StartStop & DNASeq::startTf) ) ||
        ( state == IntronF2A && (StartStop & DNASeq::startAf) ) )
  {
     return true;
  }
  return false;
}

inline bool beginSpliceStopF3(int StartStop)
{
  if  ( (StartStop & DNASeq::startTGf) || 
        (StartStop & DNASeq::startTAf) || 
        (StartStop & DNASeq::startAGf) || 
        (StartStop & DNASeq::startTCf) )
     return true;
  return false;
}

inline bool beginSpliceStopF3(int StartStop, int state)
{
  if  ( (state == IntronF3TG && (StartStop & DNASeq::startTGf) ) ||
        (state == IntronF3TA && (StartStop & DNASeq::startTAf) ) || 
        (state == IntronF3TC && (StartStop & DNASeq::startTCf) ) || 
        (state == IntronF3AG && (StartStop & DNASeq::startAGf) ) )
  {
     return true;
  }
  return false;
}

inline bool endSpliceStopF2(int StopStop, int state,  DNASeq* TheSeq)
{
  if (
      (state == IntronF2A && 
      ( (TheSeq->IsAStop("aga") && (StopStop & DNASeq::stopGAf) )||
        (TheSeq->IsAStop("agg") && (StopStop & DNASeq::stopGGf) ))  ) ||
      (state == IntronF2T && 
      ( (TheSeq->IsAStop("tga") && (StopStop & DNASeq::stopGAf) )||
        (TheSeq->IsAStop("tag") && (StopStop & DNASeq::stopAGf) )||
        (TheSeq->IsAStop("taa") && (StopStop & DNASeq::stopAAf) )||
        (TheSeq->IsAStop("tca") && (StopStop & DNASeq::stopCAf) ) ) ) 
     )
     {
       return true;
     }
     return false;
}

inline bool endSpliceStopF3(int StopStop, int state, DNASeq* TheSeq)
{
    if (
      (state == IntronF3TG && TheSeq->IsAStop("tga") && (StopStop & DNASeq::stopAf) ) ||
      (state == IntronF3TA && TheSeq->IsAStop("tag") && (StopStop & DNASeq::stopGf) ) ||
      (state == IntronF3TA && TheSeq->IsAStop("taa") && (StopStop & DNASeq::stopAf) ) ||
      (state == IntronF3TC && TheSeq->IsAStop("tca") && (StopStop & DNASeq::stopAf) ) ||
      (state == IntronF3AG && TheSeq->IsAStop("aga") && (StopStop & DNASeq::stopAf) ) ||
      (state == IntronF3AG && TheSeq->IsAStop("agg") && (StopStop & DNASeq::stopGf) ) )
  {
     return true;
  }
  return false;
}


// -------------------Reverse ---------------------------

inline bool endSpliceStopR3(int StopStop)
{
  if  ( (StopStop & DNASeq::stopAr) || (StopStop & DNASeq::stopGr))
  {
     return true;
  }
  return false;
}

inline bool endSpliceStopR3(int StopStop, int state)
{
  if  ( (state == IntronR3G && (StopStop & DNASeq::stopGr)) ||
        (state == IntronR3A && (StopStop & DNASeq::stopAr))
      )
  {
     return true;
  }
  return false;
}

inline bool endSpliceStopR2(int StopStop)
{
  if  ( (StopStop & DNASeq::stopGAr) ||  (StopStop & DNASeq::stopAGr) || 
        (StopStop & DNASeq::stopAAr) ||  (StopStop & DNASeq::stopGGr) || 
        (StopStop & DNASeq::stopCAr) )
  {
     return true;
  }
  return false;
}

inline bool endSpliceStopR2(int StopStop, int state)
{
  if  ( (state == IntronR2GA && StopStop & DNASeq::stopGAr) ||  
        (state == IntronR2AG && StopStop & DNASeq::stopAGr) || 
        (state == IntronR2AA && StopStop & DNASeq::stopAAr) || 
        (state == IntronR2CA && StopStop & DNASeq::stopCAr) || 
        (state == IntronR2GG && StopStop & DNASeq::stopGGr) )
  {
     return true;
  }
  return false;
}

inline bool beginSpliceStopR2(int StartStop, int state, DNASeq* TheSeq)
{
  if (
     (state == IntronR2AG && TheSeq->IsAStop("tag") && (StartStop & DNASeq::startTr) ) || // TAG
     (state == IntronR2GA && 
	(  ( TheSeq->IsAStop("tga") && (StartStop & DNASeq::startTr) ) ||               //TGA
	   ( TheSeq->IsAStop("aga") && (StartStop & DNASeq::startAr) ) )    )        || // AGA 
    (state == IntronR2AA  && TheSeq->IsAStop("taa") && (StartStop & DNASeq::startTr)) ||
    (state == IntronR2CA  && TheSeq->IsAStop("tca") && (StartStop & DNASeq::startTr)) ||
    (state == IntronR2GG  && TheSeq->IsAStop("agg") && (StartStop & DNASeq::startAr) ) 
    )
    {
	return true;
    }
    return false;
}


inline bool beginSpliceStopR3(int StartStop, int state, DNASeq* TheSeq)
{
  if  ( (state == IntronR3G && 
	 ( ( (StartStop & DNASeq::startTAr) && TheSeq->IsAStop("tag") ) || 
	   ( (StartStop & DNASeq::startAGr) && TheSeq->IsAStop("agg") )  )
	) ||
	(state == IntronR3A && 
	 ( ( (StartStop & DNASeq::startTGr) && TheSeq->IsAStop("tga") ) || 
	   ( (StartStop & DNASeq::startTAr) && TheSeq->IsAStop("taa") ) ||
	   ( (StartStop & DNASeq::startTCr) && TheSeq->IsAStop("tca") ) || 
	   ( (StartStop & DNASeq::startAGr) && TheSeq->IsAStop("aga") )  )
	)
      )
  {
    return true;
  }
  return false;
}

// ----------------------------------------------------------------
// A second set of macros
// ----------------------------------------------------------------
#define PICOMPEN(C,S,B,O,P) if ((C) && ISPOSSIBLE(S,B)) {               \
    BestU = PBest[Strand ? ReverseIt[O] :O]+Data.sig[DATA::S].weight[B]+(P); \
    if (isnan(maxi) || (BestU > maxi)) {maxi = BestU; best = (Strand ? ReverseIt[O] :O);}}

#define PICOMP(C,S,B,O) PICOMPEN(C,S,B,O,0.0)

#define INSERT(P)					\
  /*fprintf(stdout, "P %s, best : %s --> ", State(P).State2EGNString(), State(best).State2EGNString()); */\
  LBP[Strand ? ReverseIt[P] :P].InsertNew(best, position, maxi,PrevBP[best]);

// ----------------------------------------------------------------
// Perform DP itself: one recursive level through possible signals
// 'position' parameter is the index (start at 0)
// 'position' can be equal to the sequence length! (to manage the final signals)
// ----------------------------------------------------------------
inline void DAG::ComputeSigShifts(enum Signal::Edge Strand, DATA Data, int position)
{
  int k;
  int Data_Len = TheSeq->SeqLen;
  
  /* Note: FrameF of the first position is 0, whereas FrameR of the last position is 1 */
  int FrameF = (Strand ? ((Data_Len-position+3) % 3) : (position % 3));
  int FrameR = (Strand ? (position % 3) : ((Data_Len-position+3) % 3));

  enum Signal::Edge FStrand = (Strand ? Signal::Reverse : Signal::Forward);
  enum Signal::Edge RStrand = (Strand ? Signal::Forward : Signal::Reverse);

  double BestU, maxi = -NINFINITY;
  signed   char best = 'Z';
  bool isProkaryote = ( !strcmp(PAR.getC ("EuGene.mode"), "Prokaryote") || !strcmp(PAR.getC ("EuGene.mode"), "Prokaryote2")) ;

  // Get information on possible spliced stops
  int StopStop = TheSeq->IsStopStop(position);
  int StartStop = TheSeq->IsStartStop(position);

  Data.EstMatch = true;
  // ----------------------------------------------------------------
  // ------------------ Inits en forward ----------------------------
  // ----------------------------------------------------------------
  for (k = 0; k < 3; k++) {
    maxi = NINFINITY; best = -1;     
    
    // On commence a coder (Start), si en phase ca vient d'une UTR 5' forward
    PICOMP((FrameF == k),Start,FStrand,UTR5F);
    
    // frameshifts
    PICOMP((FrameF == (k+1)%3), Del1, FStrand, InitF1+(k+1)%3); // deletion which remove the 1st nucl of a codon
    PICOMP((FrameF == (k+2)%3), Del2, FStrand, InitF1+(k+1)%3); // deletion which remove the 2nd nucl of a codon
    PICOMP((FrameF ==  k     ), Del3, FStrand, InitF1+(k+1)%3); // deletion which remove the 3rd nucl of a codon 
    PICOMP((FrameF == (k+2)%3), Ins1, FStrand, InitF1+(k+2)%3); // insertion before the 1st nucl of a codon
    PICOMP((FrameF ==  k     ), Ins2, FStrand, InitF1+(k+2)%3); // insertion before the 2nd nucl of a codon
    PICOMP((FrameF == (k+1)%3), Ins3, FStrand, InitF1+(k+2)%3); // insertion before the 3rd nucl of a codon

    
    // On va tout droit.
      // S'il y  a un STOP en phase on ne peut continuer
    if (FrameF == k)
      LBP[Strand ? ReverseIt[InitF1+k]: InitF1+k].Update(Data.sig[DATA::Stop].weight[Signal::ForwardNo+Strand]);
    LBP[Strand ? ReverseIt[InitF1+k]: InitF1+k].Update(Data.sig[DATA::Don].weight[Signal::ForwardNo+Strand]);
    
    INSERT(InitF1+k);
  }
  
  // ----------------------------------------------------------------
  // ------------------------- Inits en reverse ---------------------
  // ----------------------------------------------------------------
  for (k = 0; k<3; k++) {
    maxi = NINFINITY; best = -1;
    
    // - on recommence a coder (Donneur) Ca vient d'un intron (no spliceable stop)
    PICOMPEN(true,Don,RStrand,IntronR1+((FrameR+3-k) % 3),
	     ((FrameR == k) ? Data.sig[DATA::Stop].weight[Signal::ReverseNo-Strand] : 0.0));
    
   if ((FrameR+3-k) % 3 == 2)
   {
      PICOMPEN(true,Don,RStrand,IntronR3G,(beginSpliceStopR3(StartStop,IntronR3G, TheSeq) ? SplicedStopPen : 0.0));
      PICOMPEN(true,Don,RStrand,IntronR3A,(beginSpliceStopR3(StartStop,IntronR3A, TheSeq) ? SplicedStopPen : 0.0));
   }
   else if ((FrameR+3-k) % 3 == 1)
   {
     PICOMPEN(true,Don,RStrand,IntronR2GA,(beginSpliceStopR2(StartStop,IntronR2GA, TheSeq) ? SplicedStopPen : 0.0));
     PICOMPEN(true,Don,RStrand,IntronR2AG,(beginSpliceStopR2(StartStop,IntronR2AG, TheSeq) ? SplicedStopPen : 0.0));
     PICOMPEN(true,Don,RStrand,IntronR2AA,(beginSpliceStopR2(StartStop,IntronR2AA, TheSeq) ? SplicedStopPen : 0.0));
     PICOMPEN(true,Don,RStrand,IntronR2CA,(beginSpliceStopR2(StartStop,IntronR2CA, TheSeq) ? SplicedStopPen : 0.0));
     PICOMPEN(true,Don,RStrand,IntronR2GG,(beginSpliceStopR2(StartStop,IntronR2GG, TheSeq) ? SplicedStopPen : 0.0));
   }
    
   
    // frameshifts
    PICOMP((FrameR ==  k     ), Del1, RStrand, InitR1+(k+2)%3); // deletion which remove the 1st nucl of a codon
    PICOMP((FrameR == (k+1)%3), Del2, RStrand, InitR1+(k+2)%3); // deletion which remove the 2nd nucl of a codon
    PICOMP((FrameR == (k+2)%3), Del3, RStrand, InitR1+(k+2)%3); // deletion which remove the 3rd nucl of a codon
    PICOMP((FrameR ==  k     ), Ins1, RStrand, InitR1+(k+1)%3); // insertion before the 1st nucl of a codon
    PICOMP((FrameR == (k+1)%3), Ins2, RStrand, InitR1+(k+1)%3); // insertion before the 2nd nucl of a codon
    PICOMP((FrameR == (k+2)%3), Ins3, RStrand, InitR1+(k+1)%3); // insertion before the 3rd nucl of a codon

    // On va tout droit.
    // S'il y  a un STOP en phase on ne peut continuer
    if ((FrameR % 3 == k)) 
      LBP[Strand ? ReverseIt[InitR1+k]: InitR1+k].Update(Data.sig[DATA::Stop].weight[Signal::ReverseNo-Strand]);
    LBP[Strand ? ReverseIt[InitR1+k]: InitR1+k].Update(Data.sig[DATA::Don].weight[Signal::ReverseNo-Strand]);
    
    INSERT(InitR1+k);
  }
  // ----------------------------------------------------------------
  // ------------------ Sngl en forward ----------------------------
  // ----------------------------------------------------------------
  for (k = 0; k < 3; k++) {
    maxi = NINFINITY; best = -1;     

    // On commence a coder (Start),si en phase ca vient d'une UTR 5' forward
    PICOMP((FrameF == k),Start,FStrand,UTR5F); 
        
    // frameshifts
    PICOMP((FrameF == (k+1)%3), Del1, FStrand, SnglF1+(k+1)%3); // deletion which remove the 1st nucl of a codon
    PICOMP((FrameF == (k+2)%3), Del2, FStrand, SnglF1+(k+1)%3); // deletion which remove the 2nd nucl of a codon
    PICOMP((FrameF ==  k     ), Del3, FStrand, SnglF1+(k+1)%3); // deletion which remove the 3rd nucl of a codon 
    PICOMP((FrameF == (k+2)%3), Ins1, FStrand, SnglF1+(k+2)%3); // insertion before the 1st nucl of a codon
    PICOMP((FrameF ==  k     ), Ins2, FStrand, SnglF1+(k+2)%3); // insertion before the 2nd nucl of a codon
    PICOMP((FrameF == (k+1)%3), Ins3, FStrand, SnglF1+(k+2)%3); // insertion before the 3rd nucl of a codon


    // ----------------------------- Prokaryote ---------------------------------------------
    // Start en phase : on peut venir d'un UIRF
    PICOMP((FrameF == k),Start, FStrand, UIRF);

    // Il ya un stop en forward, on peut venir d'une région bicodante
    // NB: pour les 6 cas suivants, inutile de tester si il y a un stop consecutif : c'est impossible, 
    PICOMP((FrameF == 1 && k == 0), Stop, FStrand, SnglF1F2); // SnglF1F2 - Stop en +2 - SnglF1
    PICOMP((FrameF == 2 && k == 0), Stop, FStrand, SnglF1F3); // SnglF1F3 - Stop en +3 - SnglF1
    PICOMP((FrameF == 0 && k == 1), Stop, FStrand, SnglF1F2); // SnglF1F2 - Stop en +1 - SnglF2
    PICOMP((FrameF == 2 && k == 1), Stop, FStrand, SnglF2F3); // SnglF2F3 - Stop en +3 - SnglF2
    PICOMP((FrameF == 0 && k == 2), Stop, FStrand, SnglF1F3); // SnglF1F3 - Stop en +1 - SnglF3
    PICOMP((FrameF == 1 && k == 2), Stop, FStrand, SnglF2F3); // SnglF2F3 - Stop en +2 - SnglF3

    // Il ya un start en reverse, on peut venir d'une région bicodante
    for (int p = 0; p < 3; p++) {
        // Si il y a un start en -p, et pas de stop en +1 si on est en phase, alors on peut venir du bicoding F1R(1+p)
        PICOMP((k == 0 && FrameR == p && !(FrameF == k && ISPOSSIBLE(Stop, FStrand)) ), Start, RStrand, SnglF1R1+p); 
        // Si il y a un start en -p, et pas de stop en +2 si on est en phase, alors on peut venir du bicoding F2R(1+p)
        PICOMP((k == 1 && FrameR == p && !(FrameF == k && ISPOSSIBLE(Stop, FStrand)) ), Start, RStrand, SnglF2R1+p); 
        // Si il y a un start en -p, et pas de stop en +3 si on est en phase, alors on peut venir du bicoding F3R(1+p)
        PICOMP((k == 2 && FrameR == p && !(FrameF == k && ISPOSSIBLE(Stop, FStrand)) ), Start, RStrand, SnglF3R1+p); 
    }

    // genes strictement adjacents
    // Note : besoin de tester si on est dans le mode procaryote. Ce n'est pas nécessaire pour les précédents cas, 
    // car les pistes UIR et Bicoding ne sont actives QUE dans le mode Procaryote.
    // Inconvenient => on fait beaucoup de tests qui s'averent inutiles dans le cas Eucaryote
    // Avantage =>  si on veut autoriser le bicoding pour le mode Eucaryote, il n'y aura pas de code a modifier. 
    if (isProkaryote)
    {
        // Un  stop en forward, suivi d'un start en forward en phase. On vient de SnglF1+k
        PICOMP((FrameF == k && ISPOSSIBLE(Stop, FStrand)), Start, FStrand, SnglF1+k);
        // un start en reverse, suivi d'un start en forward. On vient de SnglR1+FrameR
        PICOMP((FrameF == k && ISPOSSIBLE(Start, RStrand)), Start, FStrand, SnglR1+FrameR);
    }
    // -------------------------------------------------------------------------------------

    // On va tout droit.
    // S'il y  a un STOP en phase on ne peut continuer
    if ((FrameF == k))
      LBP[Strand ? ReverseIt[SnglF1+k]: SnglF1+k].Update(Data.sig[DATA::Stop].weight[Signal::ForwardNo+Strand]);
    LBP[Strand ? ReverseIt[SnglF1+k]: SnglF1+k].Update(Data.sig[DATA::Don].weight[Signal::ForwardNo+Strand]);
    
    INSERT(SnglF1+k);
  }
  // ----------------------------------------------------------------
  // ------------------------- Sngl en reverse ---------------------
  // ----------------------------------------------------------------
  for (k = 0; k<3; k++) {
    maxi = NINFINITY; best = -1;

    // On commence a coder (Stop). Ca vient d'une UTR 3' reverse
    PICOMP((FrameR == k),Stop,RStrand,UTR3R);


    // frameshifts
    PICOMP((FrameR ==  k     ), Del1, RStrand, SnglR1+(k+2)%3); // deletion which remove the 1st nucl of a codon
    PICOMP((FrameR == (k+1)%3), Del2, RStrand, SnglR1+(k+2)%3); // deletion which remove the 2nd nucl of a codon
    PICOMP((FrameR == (k+2)%3), Del3, RStrand, SnglR1+(k+2)%3); // deletion which remove the 3rd nucl of a codon
    PICOMP((FrameR ==  k     ), Ins1, RStrand, SnglR1+(k+1)%3); // insertion before the 1st nucl of a codon
    PICOMP((FrameR == (k+1)%3), Ins2, RStrand, SnglR1+(k+1)%3); // insertion before the 2nd nucl of a codon
    PICOMP((FrameR == (k+2)%3), Ins3, RStrand, SnglR1+(k+1)%3); // insertion before the 3rd nucl of a codon


    // ----------------------------- Prokaryote ---------------------------------------------
    // Stop en phase : on peut venir d'un UIRR
    PICOMP((FrameR == k),Stop,RStrand,UIRR);

    // NB: pour les 6 cas suivants, inutile de tester si il y a un stop en reverse en phase : c'est impossible, 
    // puisqu'il y a un start en reverse et que les deux etats consecutifs sont dans des frames différentes.
    PICOMP((k == 0 && FrameR == 1),Start, RStrand, SnglR1R2); // SnglR1R2, start en phase -2, SnglR1
    PICOMP((k == 0 && FrameR == 2),Start, RStrand, SnglR1R3); // SnglR1R3, start en phase -3, SnglR1
    PICOMP((k == 1 && FrameR == 0),Start, RStrand, SnglR1R2); // SnglR1R2, start en phase -1, SnglR2
    PICOMP((k == 1 && FrameR == 2),Start, RStrand, SnglR2R3); // SnglR3R2, start en phase -3, SnglR2
    PICOMP((k == 2 && FrameR == 0),Start, RStrand, SnglR1R3); // SnglR1R3, start en phase -1, SnglR3
    PICOMP((k == 2 && FrameR == 1),Start, RStrand, SnglR2R3); // SnglR2R3, start en phase -2, SnglR3

    // Si il y a un stop en +1, et pas de stop en -k si on est en phase, alors on peut venir du bicoding F1R(1+k)
    PICOMP((FrameF == 0 && !((FrameR == k) && ISPOSSIBLE(Stop, RStrand))),Stop, FStrand, SnglF1R1+k);
    // Si il y a un stop en +2, et pas de stop en -k si on est en phase, alors un peut venir du bicoding F2R(1+k)
    PICOMP((FrameF == 1 && !((FrameR == k) && ISPOSSIBLE(Stop, RStrand))),Stop, FStrand, SnglF2R1+k);
    // Si il y a un stop en +3, et pas de stop en -k si on est en phase, alors un peut venir du bicoding F3R(1+k)
    PICOMP((FrameF == 2 && !((FrameR == k) && ISPOSSIBLE(Stop, RStrand))),Stop, FStrand, SnglF3R1+k);

    // genes strictement adjacents
    // Note : besoin de tester si on est dans le mode procaryote. Ce n'est pas nécessaire pour les précédents cas, 
    // car les pistes UIR et Bicoding ne sont actives QUE dans le mode Procaryote.
    // Inconvenient => on fait beaucoup de tests qui s'averent inutiles dans le cas Eucaryote
    // Avantage =>  si on veut autoriser le bicoding pour le mode Eucaryote, il n'y aura pas de code a modifier.
    // 
    if (isProkaryote)
    {
        // un start en reverse suivi d'un stop en reverse. On peut venir de SnglR1+k
        PICOMP((FrameR == k && ISPOSSIBLE(Start, RStrand)),Stop, RStrand, SnglR1+k);
        // un stop en forward, suivi d'un stop en reverse, On peut venir de SnglF1+FrameF
        PICOMP((FrameR == k && ISPOSSIBLE(Stop, FStrand)),Stop, RStrand, SnglF1+FrameF);
    }	
    // -------------------------------------------------------------------------------------

    // On va tout droit.
    // S'il y  a un STOP en phase on ne peut continuer
    if (FrameR == k)
      LBP[Strand ? ReverseIt[SnglR1+k]: SnglR1+k].Update(Data.sig[DATA::Stop].weight[Signal::ReverseNo-Strand]);
    LBP[Strand ? ReverseIt[SnglR1+k]: SnglR1+k].Update(Data.sig[DATA::Don].weight[Signal::ReverseNo-Strand]);
    
    INSERT(SnglR1+k);
  }
  // ----------------------------------------------------------------
  // ------------------ Intrs en forward ----------------------------
  // ----------------------------------------------------------------
  for (k = 0; k < 3; k++) {
    maxi = NINFINITY; best = -1;     
    
    // On recommence a coder (Accepteur). Ca vient d'un intron (no spliceable stop)
    PICOMP(true,Acc,FStrand,IntronF1+((FrameF-k+3) % 3));

    
    if ( ((FrameF-k+3) % 3) == 1)
    {
	PICOMPEN(true,Acc,FStrand,IntronF2T, (endSpliceStopF2(StopStop,IntronF2T, TheSeq) ? SplicedStopPen : 0.0));
        PICOMPEN(true,Acc,FStrand,IntronF2A, (endSpliceStopF2(StopStop,IntronF2A, TheSeq) ? SplicedStopPen : 0.0));
    }
    else if (((FrameF-k+3) % 3) == 2)
    {
      PICOMPEN(true,Acc,FStrand,IntronF3TG, (endSpliceStopF3(StopStop,IntronF3TG, TheSeq) ? SplicedStopPen : 0.0));
      PICOMPEN(true,Acc,FStrand,IntronF3TA, (endSpliceStopF3(StopStop,IntronF3TA, TheSeq) ? SplicedStopPen : 0.0));
      PICOMPEN(true,Acc,FStrand,IntronF3TC, (endSpliceStopF3(StopStop,IntronF3TC, TheSeq) ? SplicedStopPen : 0.0));
      PICOMPEN(true,Acc,FStrand,IntronF3AG, (endSpliceStopF3(StopStop,IntronF3AG, TheSeq) ? SplicedStopPen : 0.0));
    }
    
    
    // frameshifts
    PICOMP((FrameF == (k+1)%3), Del1, FStrand, IntrF1+(k+1)%3); // deletion which remove the 1st nucl of a codon
    PICOMP((FrameF == (k+2)%3), Del2, FStrand, IntrF1+(k+1)%3); // deletion which remove the 2nd nucl of a codon
    PICOMP((FrameF ==  k     ), Del3, FStrand, IntrF1+(k+1)%3); // deletion which remove the 3rd nucl of a codon 
    PICOMP((FrameF == (k+2)%3), Ins1, FStrand, IntrF1+(k+2)%3); // insertion before the 1st nucl of a codon
    PICOMP((FrameF ==  k     ), Ins2, FStrand, IntrF1+(k+2)%3); // insertion before the 2nd nucl of a codon
    PICOMP((FrameF == (k+1)%3), Ins3, FStrand, IntrF1+(k+2)%3); // insertion before the 3rd nucl of a codon
    
    // On va tout droit.
    // S'il y  a un STOP en phase on ne peut continuer
    if (FrameF == k)
      LBP[Strand ? ReverseIt[IntrF1+k]: IntrF1+k].Update(Data.sig[DATA::Stop].weight[Signal::ForwardNo+Strand]);
    LBP[Strand ? ReverseIt[IntrF1+k]: IntrF1+k].Update(Data.sig[DATA::Don].weight[Signal::ForwardNo+Strand]);
    
    INSERT(IntrF1+k);
  }
  // ----------------------------------------------------------------
  // ------------------------- Intrs en reverse ---------------------
  // ----------------------------------------------------------------
  for (k = 0; k<3; k++) {
    maxi = NINFINITY; best = -1;
    
    // - on recommence a coder (Donneur) Ca vient d'un intron (no spliceable stop)
    PICOMPEN(true,Don,RStrand,IntronR1+((FrameR+3-k) % 3),((FrameR == k) ? Data.sig[DATA::Stop].weight[Signal::ReverseNo-Strand] : 0.0));
    
    if ((FrameR+3-k) % 3 == 2)
   {
      PICOMPEN(true,Don,RStrand,IntronR3G,(beginSpliceStopR3(StartStop,IntronR3G, TheSeq) ? SplicedStopPen : 0.0));
      PICOMPEN(true,Don,RStrand,IntronR3A,(beginSpliceStopR3(StartStop,IntronR3A, TheSeq) ? SplicedStopPen : 0.0));
   }
   else if ((FrameR+3-k) % 3 == 1)
   {
     PICOMPEN(true,Don,RStrand,IntronR2GA,(beginSpliceStopR2(StartStop,IntronR2GA, TheSeq) ? SplicedStopPen : 0.0));
     PICOMPEN(true,Don,RStrand,IntronR2AG,(beginSpliceStopR2(StartStop,IntronR2AG, TheSeq) ? SplicedStopPen : 0.0));
     PICOMPEN(true,Don,RStrand,IntronR2AA,(beginSpliceStopR2(StartStop,IntronR2AA, TheSeq) ? SplicedStopPen : 0.0));
     PICOMPEN(true,Don,RStrand,IntronR2CA,(beginSpliceStopR2(StartStop,IntronR2CA, TheSeq) ? SplicedStopPen : 0.0));
     PICOMPEN(true,Don,RStrand,IntronR2GG,(beginSpliceStopR2(StartStop,IntronR2GG, TheSeq) ? SplicedStopPen : 0.0));
   }
   

    // frameshifts
    PICOMP((FrameR ==  k     ), Del1, RStrand, IntrR1+(k+2)%3); // deletion which remove the 1st nucl of a codon
    PICOMP((FrameR == (k+1)%3), Del2, RStrand, IntrR1+(k+2)%3); // deletion which remove the 2nd nucl of a codon
    PICOMP((FrameR == (k+2)%3), Del3, RStrand, IntrR1+(k+2)%3); // deletion which remove the 3rd nucl of a codon
    PICOMP((FrameR ==  k     ), Ins1, RStrand, IntrR1+(k+1)%3); // insertion before the 1st nucl of a codon
    PICOMP((FrameR == (k+1)%3), Ins2, RStrand, IntrR1+(k+1)%3); // insertion before the 2nd nucl of a codon
    PICOMP((FrameR == (k+2)%3), Ins3, RStrand, IntrR1+(k+1)%3); // insertion before the 3rd nucl of a codon
 
    
    // On va tout droit.
    // S'il y  a un STOP en phase on ne peut continuer
    if (FrameR == k) 
      LBP[Strand ? ReverseIt[IntrR1+k]: IntrR1+k].Update(Data.sig[DATA::Stop].weight[Signal::ReverseNo-Strand]);
    LBP[Strand ? ReverseIt[IntrR1+k]: IntrR1+k].Update(Data.sig[DATA::Don].weight[Signal::ReverseNo-Strand]);
    
    INSERT(IntrR1+k);
  }
  // ----------------------------------------------------------------
  // ------------------ Terms en forward ----------------------------
  // ----------------------------------------------------------------
  for (k = 0; k < 3; k++) {
    maxi = NINFINITY; best = -1;     
    
    // On recommence a coder (Accepteur). Ca vient d'un intron
    PICOMP(true,Acc,FStrand,IntronF1+((FrameF-k+3) % 3));
    
    if ( ((FrameF-k+3) % 3) == 1)
    {
	PICOMPEN(true,Acc,FStrand,IntronF2T, (endSpliceStopF2(StopStop,IntronF2T, TheSeq) ? SplicedStopPen : 0.0));
        PICOMPEN(true,Acc,FStrand,IntronF2A, (endSpliceStopF2(StopStop,IntronF2A, TheSeq) ? SplicedStopPen : 0.0));
    }
    else if (((FrameF-k+3) % 3) == 2)
    {
      PICOMPEN(true,Acc,FStrand,IntronF3TG, (endSpliceStopF3(StopStop,IntronF3TG, TheSeq) ? SplicedStopPen : 0.0));
      PICOMPEN(true,Acc,FStrand,IntronF3TA, (endSpliceStopF3(StopStop,IntronF3TA, TheSeq) ? SplicedStopPen : 0.0));
      PICOMPEN(true,Acc,FStrand,IntronF3TC, (endSpliceStopF3(StopStop,IntronF3TC, TheSeq) ? SplicedStopPen : 0.0));
      PICOMPEN(true,Acc,FStrand,IntronF3AG, (endSpliceStopF3(StopStop,IntronF3AG, TheSeq) ? SplicedStopPen : 0.0));
    }
    
    // frameshifts
    PICOMP((FrameF == (k+1)%3), Del1, FStrand, TermF1+(k+1)%3); // deletion which remove the 1st nucl of a codon
    PICOMP((FrameF == (k+2)%3), Del2, FStrand, TermF1+(k+1)%3); // deletion which remove the 2nd nucl of a codon
    PICOMP((FrameF ==  k     ), Del3, FStrand, TermF1+(k+1)%3); // deletion which remove the 3rd nucl of a codon 
    PICOMP((FrameF == (k+2)%3), Ins1, FStrand, TermF1+(k+2)%3); // insertion before the 1st nucl of a codon
    PICOMP((FrameF ==  k     ), Ins2, FStrand, TermF1+(k+2)%3); // insertion before the 2nd nucl of a codon
    PICOMP((FrameF == (k+1)%3), Ins3, FStrand, TermF1+(k+2)%3); // insertion before the 3rd nucl of a codon
    
    // On va tout droit.
    // S'il y  a un STOP en phase on ne peut continuer
    if (FrameF == k)
      LBP[Strand ? ReverseIt[TermF1+k]: TermF1+k].Update(Data.sig[DATA::Stop].weight[Signal::ForwardNo+Strand]);
    LBP[Strand ? ReverseIt[TermF1+k]: TermF1+k].Update(Data.sig[DATA::Don].weight[Signal::ForwardNo+Strand]);
    
    INSERT(TermF1+k);
  }
  // ----------------------------------------------------------------
  // ------------------------- Terms en reverse ---------------------
  // ----------------------------------------------------------------
  for (k = 0; k<3; k++) {
    maxi = NINFINITY; best = -1;
    
    // On commence a coder (Stop). Ca vient d'une UTR 3' reverse
    PICOMP((FrameR == k),Stop,RStrand,UTR3R);
   

    // frameshifts
    PICOMP((FrameR ==  k     ), Del1, RStrand, TermR1+(k+2)%3); // deletion which remove the 1st nucl of a codon
    PICOMP((FrameR == (k+1)%3), Del2, RStrand, TermR1+(k+2)%3); // deletion which remove the 2nd nucl of a codon
    PICOMP((FrameR == (k+2)%3), Del3, RStrand, TermR1+(k+2)%3); // deletion which remove the 3rd nucl of a codon
    PICOMP((FrameR ==  k     ), Ins1, RStrand, TermR1+(k+1)%3); // insertion before the 1st nucl of a codon
    PICOMP((FrameR == (k+1)%3), Ins2, RStrand, TermR1+(k+1)%3); // insertion before the 2nd nucl of a codon
    PICOMP((FrameR == (k+2)%3), Ins3, RStrand, TermR1+(k+1)%3); // insertion before the 3rd nucl of a codon
 
    
    // On va tout droit.
    // S'il y  a un STOP en phase on ne peut continuer
    if ((FrameR == k)) 
      LBP[Strand ? ReverseIt[TermR1+k]: TermR1+k].Update(Data.sig[DATA::Stop].weight[Signal::ReverseNo-Strand]);
    LBP[Strand ? ReverseIt[TermR1+k]: TermR1+k].Update(Data.sig[DATA::Don].weight[Signal::ReverseNo-Strand]);
    
    INSERT(TermR1+k);
  }
  // ----------------------------------------------------------------
  // ------------------------ Intergenique --------------------------
  // ----------------------------------------------------------------
  // Ca peut venir d'une fin de 3' direct ou de 5' reverse
  maxi = NINFINITY; best = -1;
  
  // From 5' reverse
  PICOMP(true,tStartNpc,RStrand, RnaR);
  PICOMP(true,tStart,RStrand, UTR5R);

  // From 3' direct
  PICOMP(true,tStopNpc, FStrand, RnaF);
  PICOMP(true,tStop,    FStrand, UTR3F);

  // On reste intergenique
  // et les transstartNO/TransstopNO ???
  
  INSERT(InterGen);
  // ----------------------------------------------------------------
  // ---------------------- UTR 5' direct ---------------------------
  // ----------------------------------------------------------------
  maxi = NINFINITY; best = -1;
  
  // On vient de l'intergenique. 
  PICOMP(true,tStart,FStrand, InterGen);
  // On peut venir aussi d'un intron d'UTR.
  PICOMP((Data.EstMatch),Acc,FStrand, IntronU5F);
  
  // On reste 5' direct. On ne prend pas le Start eventuel.
  LBP[Strand ? ReverseIt[UTR5F]: UTR5F].Update(Data.sig[DATA::Start].weight[Signal::ForwardNo+Strand]);
  
  INSERT(UTR5F);

  // ----------------------------------------------------------------
  // ---------------------- Rna Forward ---------------------------
  // ----------------------------------------------------------------
  maxi = NINFINITY; best = -1;

  // From Intergenic
  PICOMP(true, tStartNpc, FStrand, InterGen);
  // From UIRF
  PICOMP(true, tStartNpc, FStrand, UIRF);

  //LBP[Strand ? ReverseIt[RnaF]: RnaF].Update(Data.sig[DATA::Start].weight[Signal::ForwardNo+Strand]);

  INSERT(RnaF);

  

  // ----------------------------------------------------------------
  // ---------------------- Rna Reverse ---------------------------
  // ----------------------------------------------------------------
  maxi = NINFINITY; best = -1;

  // From Intergenic 
  PICOMP(true,tStopNpc, RStrand, InterGen);
  // From UIR Reverse
  PICOMP(true,tStopNpc, RStrand, UIRR);

  //LBP[Strand ? ReverseIt[RnaR]: RnaR].Update(Data.sig[DATA::Start].weight[Signal::ReverseNo-Strand]);
  INSERT(RnaR);

  // ----------------------------------------------------------------
  // ------------------- Introns d'UTR5F ----------------------------
  // ----------------------------------------------------------------
  maxi = NINFINITY; best = -1;
  
  // on quitte l'UTR 
  PICOMP((Data.EstMatch),Don,FStrand, UTR5F);
  
  // On reste intronique
  LBP[Strand ? ReverseIt[IntronU5F]: IntronU5F].Update(Data.sig[DATA::Acc].weight[Signal::ForwardNo+Strand]);
  
  INSERT(IntronU5F);
  // ----------------------------------------------------------------
  // ------------------- Introns d'UTR3F ----------------------------
  // ----------------------------------------------------------------
  maxi = NINFINITY; best = -1;
  
  // on quitte l'UTR 
  PICOMP((Data.EstMatch),Don,FStrand, UTR3F);
  
  // On reste intronique
  LBP[Strand ? ReverseIt[IntronU3F]: IntronU3F].Update(Data.sig[DATA::Acc].weight[Signal::ForwardNo+Strand]);
  
  INSERT(IntronU3F);
  // ----------------------------------------------------------------
  // ---------------------- UTR 3' direct ---------------------------
  // ----------------------------------------------------------------
  maxi = NINFINITY; best = -1;
  
  // Ca vient d'un Term ou d'un Sngl direct + STOP
  PICOMP(true,Stop,FStrand, TermF1+FrameF);
  PICOMP(true,Stop,FStrand, SnglF1+FrameF);
  // On peut venir aussi d'un intron d'UTR.
  PICOMP((Data.EstMatch),Acc,FStrand, IntronU3F);
  
  // On reste 3' direct
  
  INSERT(UTR3F);
  // ----------------------------------------------------------------
  // ----------------------- UTR 5'reverse --------------------------
  // ----------------------------------------------------------------
  maxi = NINFINITY; best = -1;
  
  // Ca vient d'un Sngl ou Init reverse + START
  PICOMP(true,Start,RStrand, InitR1+FrameR);
  PICOMP(true,Start,RStrand, SnglR1+FrameR);
  // On peut venir aussi d'un intron d'UTR.
  PICOMP((Data.EstMatch),Don,RStrand, IntronU5R);
  
  // On reste 5' reverse
  LBP[Strand ? ReverseIt[UTR5R]: UTR5R].Update(Data.sig[DATA::Start].weight[Signal::ReverseNo-Strand]);
  
  INSERT(UTR5R);
  // ----------------------------------------------------------------
  // ------------------- Introns d'UTR5R ----------------------------
  // ----------------------------------------------------------------
  maxi = NINFINITY; best = -1;
  
  // On quitte une UTR5R
  PICOMP((Data.EstMatch),Acc,RStrand, UTR5R);
  
  // On reste intronique
  LBP[Strand ? ReverseIt[IntronU5R]: IntronU5R].Update(Data.sig[DATA::Acc].weight[Signal::ReverseNo-Strand]);
  
  INSERT(IntronU5R);
  // ----------------------------------------------------------------
  // ------------------- Introns d'UTR3R ----------------------------
  // ----------------------------------------------------------------
  maxi = NINFINITY; best = -1;
  
  // On quitte une UTR5R
  PICOMP((Data.EstMatch),Acc,RStrand, UTR3R);
  
  // On reste intronique
  LBP[Strand ? ReverseIt[IntronU3R]: IntronU3R].Update(Data.sig[DATA::Acc].weight[Signal::ReverseNo-Strand]);
  
  INSERT(IntronU3R);
  // ----------------------------------------------------------------
  // ----------------------- UTR 3'reverse --------------------------
  // ----------------------------------------------------------------
  maxi = NINFINITY; best = -1;
  
  // On demarre depuis l'intergenique
  PICOMP(true,tStop,RStrand, InterGen);
  // On peut venir aussi d'un intron d'UTR.
  PICOMP((Data.EstMatch),Don,RStrand, IntronU3R);
  
  // On reste 3' reverse
  
  INSERT(UTR3R);
  // ----------------------------------------------------------------
  // ---------------- Introns de phase k forward --------------------
  // ----------------------------------------------------------------
  for (k = 0; k<3; k++) {
    maxi = NINFINITY; best = -1;
    
    // - on quitte un Init ou un Intr
    // no spliceable stop: 
    if (! ( ( k == 1 && beginSpliceStopF2(StartStop) ) || 
	    ( k == 2 && beginSpliceStopF3(StartStop) ) )   )  {
      PICOMPEN(true,Don,FStrand, InitF1+((FrameF-k+3) % 3),
	       (k == 0 ? Data.sig[DATA::Stop].weight[Signal::ForwardNo+Strand] : 0.0));
      PICOMPEN(true,Don,FStrand, IntrF1+((FrameF-k+3) % 3),
	       (k == 0 ? Data.sig[DATA::Stop].weight[Signal::ForwardNo+Strand] : 0.0));
    }
    
    // On reste intronique
    LBP[Strand ? ReverseIt[IntronF1+k]: IntronF1+k].Update(Data.sig[DATA::Acc].weight[Signal::ForwardNo+Strand]);
    
    INSERT(IntronF1+k);
  }
  // ----------------------------------------------------------------
  // ----------------- Introns forward speciaux ---------------------
  // ----------------------------------------------------------------
  
  maxi = NINFINITY; best = -1;
  k = 1;
  if (beginSpliceStopF2(StartStop, IntronF2T)) {
    // - on quitte un Init ou un Intr
    PICOMP(true,Don,FStrand, InitF1+((FrameF-k+3) % 3));
    PICOMP(true,Don,FStrand, IntrF1+((FrameF-k+3) % 3));
  }
  // On reste intronique
  LBP[Strand ? ReverseIt[IntronF2T]: IntronF2T].Update(Data.sig[DATA::Acc].weight[Signal::ForwardNo+Strand]);
  INSERT(IntronF2T);
  
  maxi = NINFINITY; best = -1;
  k = 1;
  if (beginSpliceStopF2(StartStop, IntronF2A)) {
    // - on quitte un Init ou un Intr
    PICOMP(true,Don,FStrand, InitF1+((FrameF-k+3) % 3));
    PICOMP(true,Don,FStrand, IntrF1+((FrameF-k+3) % 3));
  }
  // On reste intronique
  LBP[Strand ? ReverseIt[IntronF2A]: IntronF2A].Update(Data.sig[DATA::Acc].weight[Signal::ForwardNo+Strand]);
  INSERT(IntronF2A);
  
  maxi = NINFINITY; best = -1;
  k = 2;
  if (beginSpliceStopF3(StartStop,IntronF3TG)) {
    // - on quitte un Init ou un Intr
    PICOMP(true,Don,FStrand, InitF1+((FrameF-k+3) % 3));
    PICOMP(true,Don,FStrand, IntrF1+((FrameF-k+3) % 3));
  }
  // On reste intronique
  LBP[Strand ? ReverseIt[IntronF3TG]: IntronF3TG].Update(Data.sig[DATA::Acc].weight[Signal::ForwardNo+Strand]);
  INSERT(IntronF3TG);
  
  
  maxi = NINFINITY; best = -1;
  k = 2;
  if (beginSpliceStopF3(StartStop,IntronF3TA)) {
    // - on quitte un Init ou un Intr
    PICOMP(true,Don,FStrand, InitF1+((FrameF-k+3) % 3));
    PICOMP(true,Don,FStrand, IntrF1+((FrameF-k+3) % 3));
  }
  // On reste intronique
  LBP[Strand ? ReverseIt[IntronF3TA]: IntronF3TA].Update(Data.sig[DATA::Acc].weight[Signal::ForwardNo+Strand]);
  INSERT(IntronF3TA);
  
  maxi = NINFINITY; best = -1;
  k = 2;
  if (beginSpliceStopF3(StartStop,IntronF3TC)) {
    // - on quitte un Init ou un Intr
    PICOMP(true,Don,FStrand, InitF1+((FrameF-k+3) % 3));
    PICOMP(true,Don,FStrand, IntrF1+((FrameF-k+3) % 3));
  }
  // On reste intronique
  LBP[Strand ? ReverseIt[IntronF3TC]: IntronF3TC].Update(Data.sig[DATA::Acc].weight[Signal::ForwardNo+Strand]);
  INSERT(IntronF3TC);
  
  maxi = NINFINITY; best = -1;
  k = 2;
  if (beginSpliceStopF3(StartStop,IntronF3AG)) {
    // - on quitte un Init ou un Intr
    PICOMP(true,Don,FStrand, InitF1+((FrameF-k+3) % 3));
    PICOMP(true,Don,FStrand, IntrF1+((FrameF-k+3) % 3));
  }
  // On reste intronique
  LBP[Strand ? ReverseIt[IntronF3AG]: IntronF3AG].Update(Data.sig[DATA::Acc].weight[Signal::ForwardNo+Strand]);
  INSERT(IntronF3AG);
  
 
  // ----------------------------------------------------------------
  // ----------------- Introns de phase -k reverse ------------------
  // ----------------------------------------------------------------
  for (k = 0; k<3; k++) {
    maxi = NINFINITY; best = -1;
    
    // On quitte un Intr ou un Term reverse
    // no spliceable stop: 
    if (! ( ( (k == 2) && endSpliceStopR3(StopStop) ) || 
	    ( (k == 1) && endSpliceStopR2(StopStop) ) ) ) 
    {
      PICOMP(true,Acc,RStrand, IntrR1+((FrameR+3-k) % 3));
      PICOMP(true,Acc,RStrand, TermR1+((FrameR+3-k) % 3));
    }
    
    // On reste intronique
    LBP[Strand ? ReverseIt[IntronR1+k]: IntronR1+k].Update(Data.sig[DATA::Acc].weight[Signal::ReverseNo-Strand]);
    
    INSERT(IntronR1+k);
  }
  // ----------------------------------------------------------------
  // ----------------- Introns reverse speciaux ---------------------
  // ----------------------------------------------------------------
  
  maxi = NINFINITY; best = -1;
  k = 2;
  if (endSpliceStopR3(StopStop,IntronR3G)) {
    // - on quitte un Intr ou un Term
    PICOMP(true,Acc,RStrand, IntrR1+((FrameR+3-k) % 3));
    PICOMP(true,Acc,RStrand, TermR1+((FrameR+3-k) % 3));
  }
  // On reste intronique
  LBP[Strand ? ReverseIt[IntronR3G]: IntronR3G].Update(Data.sig[DATA::Acc].weight[Signal::ReverseNo-Strand]);
  INSERT(IntronR3G);
  
  maxi = NINFINITY; best = -1;
  k = 2;
  if (endSpliceStopR3(StopStop,IntronR3A)) {
    // - on quitte un Intr ou un Term
    PICOMP(true,Acc,RStrand, IntrR1+((FrameR+3-k) % 3));
    PICOMP(true,Acc,RStrand, TermR1+((FrameR+3-k) % 3));
  }
  // On reste intronique
  LBP[Strand ? ReverseIt[IntronR3A]: IntronR3A].Update(Data.sig[DATA::Acc].weight[Signal::ReverseNo-Strand]);
  INSERT(IntronR3A);
  
 
  maxi = NINFINITY; best = -1;
  k = 1;
  if (endSpliceStopR2(StopStop,IntronR2GA)) {
    // - on quitte un Intr ou un Term
    PICOMP(true,Acc,RStrand, IntrR1+((FrameR+3-k) % 3));
    PICOMP(true,Acc,RStrand, TermR1+((FrameR+3-k) % 3));
  }
  // On reste intronique
  LBP[Strand ? ReverseIt[IntronR2GA]: IntronR2GA].Update(Data.sig[DATA::Acc].weight[Signal::ReverseNo-Strand]);
  INSERT(IntronR2GA); 
  
  maxi = NINFINITY; best = -1;
  k = 1;
  if (endSpliceStopR2(StopStop,IntronR2AG)) {
    // - on quitte un Intr ou un Term
    PICOMP(true,Acc,RStrand, IntrR1+((FrameR+3-k) % 3));
    PICOMP(true,Acc,RStrand, TermR1+((FrameR+3-k) % 3));
  }
  // On reste intronique
  LBP[Strand ? ReverseIt[IntronR2AG]: IntronR2AG].Update(Data.sig[DATA::Acc].weight[Signal::ReverseNo-Strand]);
  INSERT(IntronR2AG); 
  
  maxi = NINFINITY; best = -1;
  k = 1;
  if (endSpliceStopR2(StopStop,IntronR2AA)) {
    // - on quitte un Intr ou un Term
    PICOMP(true,Acc,RStrand, IntrR1+((FrameR+3-k) % 3));
    PICOMP(true,Acc,RStrand, TermR1+((FrameR+3-k) % 3));
  }
  // On reste intronique
  LBP[Strand ? ReverseIt[IntronR2AA]: IntronR2AA].Update(Data.sig[DATA::Acc].weight[Signal::ReverseNo-Strand]);
  INSERT(IntronR2AA); 
  
  maxi = NINFINITY; best = -1;
  k = 1;
  if (endSpliceStopR2(StopStop,IntronR2CA)) {
    // - on quitte un Intr ou un Term
    PICOMP(true,Acc,RStrand, IntrR1+((FrameR+3-k) % 3));
    PICOMP(true,Acc,RStrand, TermR1+((FrameR+3-k) % 3));
  }
  // On reste intronique
  LBP[Strand ? ReverseIt[IntronR2CA]: IntronR2CA].Update(Data.sig[DATA::Acc].weight[Signal::ReverseNo-Strand]);
  INSERT(IntronR2CA);
  
  maxi = NINFINITY; best = -1;
  k = 1;
  if (endSpliceStopR2(StopStop,IntronR2GG)) {
    // - on quitte un Intr ou un Term
    PICOMP(true,Acc,RStrand, IntrR1+((FrameR+3-k) % 3));
    PICOMP(true,Acc,RStrand, TermR1+((FrameR+3-k) % 3));
  }
  // On reste intronique
  LBP[Strand ? ReverseIt[IntronR2GG]: IntronR2GG].Update(Data.sig[DATA::Acc].weight[Signal::ReverseNo-Strand]);
  INSERT(IntronR2GG); 
  

  // ----------------------------------------------------------------
  // ------------------ SnglF1F2 ------------------------------------
  // ----------------------------------------------------------------
  maxi = NINFINITY; best = -1;
  PICOMP( (FrameF == 1), Start, FStrand, SnglF1); // je suis en phase 2, je viens de phase1
  PICOMP( (FrameF == 0), Start, FStrand, SnglF2); // je suis en phase 1, je viens de phase2
  // FrameShift??

    // On va tout droit.
	// S'il y  a un STOP en phase on ne peut continuer
    if (FrameF == 0 || FrameF == 1)
      LBP[Strand ? ReverseIt[SnglF1F2]: SnglF1F2].Update(Data.sig[DATA::Stop].weight[Signal::ForwardNo+Strand]);

  INSERT(SnglF1F2);
  // ----------------------------------------------------------------
  // ------------------ SnglF1F3 ------------------------------------
  // ----------------------------------------------------------------
  maxi = NINFINITY; best = -1;
  PICOMP( (FrameF == 2), Start, FStrand, SnglF1); // je suis en phase 3, je viens de phase1
  PICOMP( (FrameF == 0), Start, FStrand, SnglF3); // je suis en phase 1, je viens de phase2
  // FrameShift??

    // On va tout droit.
  // S'il y  a un STOP en phase on ne peut continuer
  if (FrameF == 0 || FrameF == 2)
      LBP[Strand ? ReverseIt[SnglF1F3]: SnglF1F3].Update(Data.sig[DATA::Stop].weight[Signal::ForwardNo+Strand]);

  INSERT(SnglF1F3);
  // ----------------------------------------------------------------
  // ------------------ SnglF2F3 ------- ----------------------------
  // ----------------------------------------------------------------
  maxi = NINFINITY; best = -1;
  PICOMP( (FrameF == 2), Start, FStrand, SnglF2); // je suis en phase 3, je viens de phase2
  PICOMP( (FrameF == 1), Start, FStrand, SnglF3); // je suis en phase 2, je viens de phase3
  // FrameShift??

  // On va tout droit.
  // S'il y  a un STOP en phase on ne peut continuer
  if (FrameF == 1 || FrameF == 2)
      LBP[Strand ? ReverseIt[SnglF2F3]: SnglF2F3].Update(Data.sig[DATA::Stop].weight[Signal::ForwardNo+Strand]);

  INSERT(SnglF2F3);
  // ----------------------------------------------------------------
  // ------------------ SnglR1R2 ------------------------------------
  // ----------------------------------------------------------------
  maxi = NINFINITY; best = -1;
  PICOMP( (FrameR == 0), Stop, RStrand, SnglR2); // je suis en phase 1, je viens de phase2
  PICOMP( (FrameR == 1), Stop, RStrand, SnglR1); // je suis en phase 2, je viens de phase1
  // FrameShift??

  // On va tout droit.
  // S'il y  a un STOP en phase on ne peut continuer
  if (FrameR == 0 || FrameR == 1)
	LBP[Strand ? ReverseIt[SnglR1R2]: SnglR1R2].Update(Data.sig[DATA::Stop].weight[Signal::ReverseNo-Strand]);

  INSERT(SnglR1R2);
  // ----------------------------------------------------------------
  // ------------------ SnglR1R3 ------------------------------------
  // ----------------------------------------------------------------
  maxi = NINFINITY; best = -1;
  PICOMP( (FrameR == 0), Stop, RStrand, SnglR3); // je suis en phase 1, je viens de phase3
  PICOMP( (FrameR == 2), Stop, RStrand, SnglR1); // je suis en phase 3, je viens de phase1
  // FrameShift??

  // On va tout droit.
  // S'il y  a un STOP en phase on ne peut continuer
   if (FrameR == 0 || FrameR == 2)
	LBP[Strand ? ReverseIt[SnglR1R3]: SnglR1R3].Update(Data.sig[DATA::Stop].weight[Signal::ReverseNo-Strand]);

  INSERT(SnglR1R3);
  // ----------------------------------------------------------------
  // ------------------ SnglR2R3 ------------------------------------
  // ----------------------------------------------------------------
  maxi = NINFINITY; best = -1;
  PICOMP( (FrameR == 1), Stop, RStrand, SnglR3); // je suis en phase 2, je viens de phase3
  PICOMP( (FrameR == 2), Stop, RStrand, SnglR2); // je suis en phase 3, je viens de phase2
  // FrameShift??

  // On va tout droit.
  // S'il y  a un STOP en phase on ne peut continuer
  if (FrameR == 1 || FrameR == 2)
	LBP[Strand ? ReverseIt[SnglR2R3]: SnglR2R3].Update(Data.sig[DATA::Stop].weight[Signal::ReverseNo-Strand]);

  INSERT(SnglR2R3);

  // ----------------------------------------------------------------
  // ------------------ SnglF1R1  -----------------------------------
  // ----------------------------------------------------------------
   maxi = NINFINITY; best = -1; 
   PICOMP(( (!( (FrameF == 0) && ISPOSSIBLE(Stop, FStrand))) && (FrameR == 0) ), Stop,  RStrand, SnglF1);
   PICOMP(( (!( (FrameR == 0) && ISPOSSIBLE(Stop, RStrand))) && (FrameF == 0) ), Start, FStrand, SnglR1);
   // FrameShift??

    // On va tout droit.
  // S'il y  a un STOP en phase on ne peut continuer
  if (FrameF == 0)
    LBP[Strand ? ReverseIt[SnglF1R1]: SnglF1R1].Update(Data.sig[DATA::Stop].weight[Signal::ForwardNo+Strand]);
  if (FrameR == 0)
	LBP[Strand ? ReverseIt[SnglF1R1]: SnglF1R1].Update(Data.sig[DATA::Stop].weight[Signal::ReverseNo-Strand]);

  INSERT(SnglF1R1);

  // ----------------------------------------------------------------
  // ------------------ SnglF1R2  -----------------------------------
  // ----------------------------------------------------------------
   maxi = NINFINITY; best = -1;
   PICOMP(( (!( (FrameF == 0) && ISPOSSIBLE(Stop, FStrand))) && FrameR == 1), Stop,  RStrand, SnglF1);
   PICOMP(( (!( (FrameR == 1) && ISPOSSIBLE(Stop, RStrand))) && FrameF == 0), Start, FStrand, SnglR2);
   // FrameShift??

  // On va tout droit.
  // S'il y  a un STOP en phase on ne peut continuer
  if (FrameF == 0)
    LBP[Strand ? ReverseIt[SnglF1R2]: SnglF1R2].Update(Data.sig[DATA::Stop].weight[Signal::ForwardNo+Strand]);
  if (FrameR == 1)
	LBP[Strand ? ReverseIt[SnglF1R2]: SnglF1R2].Update(Data.sig[DATA::Stop].weight[Signal::ReverseNo-Strand]);

  INSERT(SnglF1R2);

  // ----------------------------------------------------------------
  // ------------------ SnglF1R3  -----------------------------------
  // ----------------------------------------------------------------
   maxi = NINFINITY; best = -1;
   PICOMP(( (!( (FrameF == 0) && ISPOSSIBLE(Stop, FStrand))) && FrameR == 2), Stop,  RStrand, SnglF1);
   PICOMP(( (!( (FrameR == 2) && ISPOSSIBLE(Stop, RStrand))) && FrameF == 0), Start, FStrand, SnglR3);
   // FrameShift??

  // On va tout droit.
  // S'il y  a un STOP en phase on ne peut continuer
  if (FrameF == 0)
    LBP[Strand ? ReverseIt[SnglF1R3]: SnglF1R3].Update(Data.sig[DATA::Stop].weight[Signal::ForwardNo+Strand]);
  if (FrameR == 2)
	LBP[Strand ? ReverseIt[SnglF1R3]: SnglF1R3].Update(Data.sig[DATA::Stop].weight[Signal::ReverseNo-Strand]);

  INSERT(SnglF1R3);


  // ----------------------------------------------------------------
  // ------------------ SnglF2R1  -----------------------------------
  // ----------------------------------------------------------------
   maxi = NINFINITY; best = -1;
   PICOMP(((!( (FrameF == 1) && ISPOSSIBLE(Stop, FStrand))) && FrameR == 0), Stop,  RStrand, SnglF2);
   PICOMP(((!( (FrameR == 0) && ISPOSSIBLE(Stop, RStrand))) && FrameF == 1), Start, FStrand, SnglR1);
   // FrameShift??

  // On va tout droit.
  // S'il y  a un STOP en phase on ne peut continuer
  if (FrameF == 1)
    LBP[Strand ? ReverseIt[SnglF2R1]: SnglF2R1].Update(Data.sig[DATA::Stop].weight[Signal::ForwardNo+Strand]);
  if (FrameR == 0)
	LBP[Strand ? ReverseIt[SnglF2R1]: SnglF2R1].Update(Data.sig[DATA::Stop].weight[Signal::ReverseNo-Strand]);

  INSERT(SnglF2R1);

  // ----------------------------------------------------------------
  // ------------------ SnglF2R2  -----------------------------------
  // ----------------------------------------------------------------
   maxi = NINFINITY; best = -1;
   PICOMP(( (!( (FrameF == 1) && ISPOSSIBLE(Stop, FStrand))) && FrameR == 1), Stop,  RStrand, SnglF2);
   PICOMP(( (!( (FrameR == 1) && ISPOSSIBLE(Stop, RStrand))) && FrameF == 1), Start, FStrand, SnglR2);
   // FrameShift??
  // On va tout droit.
  // S'il y  a un STOP en phase on ne peut continuer
  if (FrameF == 1)
    LBP[Strand ? ReverseIt[SnglF2R2]: SnglF2R2].Update(Data.sig[DATA::Stop].weight[Signal::ForwardNo+Strand]);
  if (FrameR == 1)
	LBP[Strand ? ReverseIt[SnglF2R2]: SnglF2R2].Update(Data.sig[DATA::Stop].weight[Signal::ReverseNo-Strand]);

  INSERT(SnglF2R2);

  // ----------------------------------------------------------------
  // ------------------ SnglF2R3  -----------------------------------
  // ----------------------------------------------------------------
   maxi = NINFINITY; best = -1;
   PICOMP(((!( (FrameF == 1) && ISPOSSIBLE(Stop, FStrand))) && FrameR == 2), Stop,  RStrand, SnglF2);
   PICOMP(((!( (FrameR == 2) && ISPOSSIBLE(Stop, RStrand))) && FrameF == 1), Start, FStrand, SnglR3);
   // FrameShift??

  // On va tout droit.
  // S'il y  a un STOP en phase on ne peut continuer
    if (FrameF == 1)
    {
    LBP[Strand ? ReverseIt[SnglF2R3]: SnglF2R3].Update(Data.sig[DATA::Stop].weight[Signal::ForwardNo+Strand]);
    }
    if (FrameR == 2)
		LBP[Strand ? ReverseIt[SnglF2R3]: SnglF2R3].Update(Data.sig[DATA::Stop].weight[Signal::ReverseNo-Strand]);

  INSERT(SnglF2R3);


  // ----------------------------------------------------------------
  // ------------------ SnglF3R1  -----------------------------------
  // ----------------------------------------------------------------
   maxi = NINFINITY; best = -1;
   PICOMP( ((!( (FrameF == 2) && ISPOSSIBLE(Stop, FStrand))) && FrameR == 0), Stop,  RStrand, SnglF3);
   PICOMP( ((!( (FrameR == 0) && ISPOSSIBLE(Stop, RStrand))) && FrameF == 2), Start, FStrand, SnglR1);
   // FrameShift??
  // On va tout droit.

    // S'il y  a un STOP en phase on ne peut continuer
  if (FrameF == 2)
    LBP[Strand ? ReverseIt[SnglF3R1]: SnglF3R1].Update(Data.sig[DATA::Stop].weight[Signal::ForwardNo+Strand]);
  if (FrameR == 0)
	LBP[Strand ? ReverseIt[SnglF3R1]: SnglF3R1].Update(Data.sig[DATA::Stop].weight[Signal::ReverseNo-Strand]);

  INSERT(SnglF3R1);

  // ----------------------------------------------------------------
  // ------------------ SnglF3R2  -----------------------------------
  // ----------------------------------------------------------------
   maxi = NINFINITY; best = -1;
   PICOMP(((!( (FrameF == 2) && ISPOSSIBLE(Stop, FStrand))) && FrameR == 1), Stop,  RStrand, SnglF3);
   PICOMP(((!( (FrameR == 1) && ISPOSSIBLE(Stop, RStrand))) && FrameF == 2), Start, FStrand, SnglR2);
   // FrameShift??

  // On va tout droit.
  // S'il y  a un STOP en phase on ne peut continuer
  if (FrameF == 2)
    LBP[Strand ? ReverseIt[SnglF3R2]: SnglF3R2].Update(Data.sig[DATA::Stop].weight[Signal::ForwardNo+Strand]);
  if (FrameR == 1)
	LBP[Strand ? ReverseIt[SnglF3R2]: SnglF3R2].Update(Data.sig[DATA::Stop].weight[Signal::ReverseNo-Strand]);

  INSERT(SnglF3R2);


  // ----------------------------------------------------------------
  // ------------------ SnglF3R3  -----------------------------------
  // ----------------------------------------------------------------
   maxi = NINFINITY; best = -1;
   PICOMP(((!( (FrameF == 2) && ISPOSSIBLE(Stop, FStrand))) && FrameR == 2), Stop,  RStrand, SnglF3);
   PICOMP(((!( (FrameR == 2) && ISPOSSIBLE(Stop, RStrand))) && FrameF == 2), Start, FStrand, SnglR3);
   // FrameShift??

  // On va tout droit.
  // S'il y  a un STOP en phase on ne peut continuer
   if (FrameF == 2)
    LBP[Strand ? ReverseIt[SnglF3R3]: SnglF3R3].Update(Data.sig[DATA::Stop].weight[Signal::ForwardNo+Strand]);
  if (FrameR == 2)
	LBP[Strand ? ReverseIt[SnglF3R3]: SnglF3R3].Update(Data.sig[DATA::Stop].weight[Signal::ReverseNo-Strand]);


  INSERT(SnglF3R3);

  // ----------------------------------------------------------------
  // ---------------------- UIR forward ---------------------------
  // ----------------------------------------------------------------
  maxi = NINFINITY; best = -1;
  
  // From a Sngl
  PICOMP(true,Stop,FStrand, SnglF1+FrameF);
  // From a RnaF
  PICOMP(true,tStopNpc, FStrand, RnaF);

  INSERT(UIRF);

  // ----------------------------------------------------------------
  // ---------------------- UIR reverse ---------------------------
  // ----------------------------------------------------------------
  maxi = NINFINITY; best = -1;
  
  // Ca vient d'un Sngl + START
  PICOMP(true, Start, RStrand, SnglR1+FrameR);
  // From a RNAR
  PICOMP(true, tStartNpc, RStrand, RnaR);

  INSERT(UIRR);
  
}

// ----------------------------------------------------------------
//  Shortest Path Algorithm with length constraint
// ----------------------------------------------------------------
void DAG :: ShortestPathAlgoForward (int position, DATA Data)
{
  //fprintf(stderr, "POSITION: %d\n", position);
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
