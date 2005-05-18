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
// File:     BackP.cc
// Contents: Definitions for a linear time shortest path with constraint alg.
// ------------------------------------------------------------------

#include <stdio.h>
#include <math.h>

#include "BackP.h"

#include "System.h"

// ----------------------------------------------------------------
// Default constructor.
// ----------------------------------------------------------------
BackPoint :: BackPoint  ()
{
  State = -1;
  StartPos = INITIALSHIFT;
  Cost = 0.0;
  Additional = 0.0;
  Next = Prev = this;
  Origin = NULL;
  Optimal = false;
}
// ----------------------------------------------------------------
//  Default constructor.
// ----------------------------------------------------------------
BackPoint :: BackPoint  (char state, int pos, double cost)
{
  State = state;
  StartPos = pos;
  Cost = cost;
  Additional = 0.0;
  Next = Prev = Origin = NULL;
  Optimal = false;
}

// ----------------------------------------------------------------
//  Default destructor.
// ----------------------------------------------------------------
BackPoint :: ~BackPoint  ()
{
  Prev->Next = Next;
  Next->Prev = Prev;
}

// ----------------------------------------------------------------
// Prints the BackPoint contents
// ----------------------------------------------------------------
void BackPoint :: Print()
{
  printf("pos = %d, state = %d\n", StartPos,State);
}
// ----------------------------------------------------------------
// Track creator
// ----------------------------------------------------------------
Track :: Track()
{
  NumBP = 0;
  Optimal = 0.0;
  OptPos =0;
}
// ----------------------------------------------------------------
// Track destructor
// ----------------------------------------------------------------
Track :: ~Track()
{
  Zap();
}
// ----------------------------------------------------------------
// Zap  the path data of a  whole track
// ----------------------------------------------------------------
void Track :: Zap()
{
  BackPoint *Dead = Path.Next,*Skip;

  while (Dead != &Path) {
    Skip = Dead->Next;
    delete Dead;
    Dead = Skip;
  }
}
// ----------------------------------------------------------------
// Insert  a new backpoint
// ----------------------------------------------------------------
void Track :: InsertNew(char state, int pos, double cost, BackPoint *Or)
{
  BackPoint *It = Path.Next;
  
  //  if (cost > NINFINITY) {
  //  if (cost > Optimal-PenD.MaxDelta) {
  //  if (cost > Optimal-PenD.GetDelta(pos-OptPos)) {
  //  if (cost > Path.Next->Cost+Path.Next->Additional-PenD.GetDelta(pos-Path.Next->StartPos)) {
  if (cost > Optimal-PenD.GetDelta(pos-OptPos) && 
      cost > Path.Next->Cost+Path.Next->Additional-PenD.GetDelta(pos-Path.Next->StartPos)) {
    NumBP++;
    It =  new BackPoint(state,pos,cost);
    if (cost > Optimal) { 
      Optimal =cost; 
      It->Optimal = true;
      OptPos = pos;
    }
    It->Next = Path.Next;
    It->Prev = Path.Next->Prev;
    It->Origin = Or;
    Path.Next->Prev = It;
    Path.Next = It;
  }
}
// ----------------------------------------------------------------
// Insert  a new backpoint
// ----------------------------------------------------------------
void Track :: ForceNew(char state, int pos, double cost, BackPoint *Or)
{
  BackPoint *It = Path.Next;
  
  NumBP++;
  It =  new BackPoint(state,pos,cost);
  if (cost > Optimal) { Optimal =cost; It->Optimal = true;}
  It->Next = Path.Next;
  It->Prev = Path.Next->Prev;
  It->Origin = Or;
  Path.Next->Prev = It;
  Path.Next = It;
}
// ----------------------------------------------------------------
// Returns the best BackPoint and the BackPoint is at least len
// nuc. far from pos
// ----------------------------------------------------------------
BackPoint *Track :: BestUsable(int pos, double *cost, int Forward, int pen)
{
  BackPoint *BestBP = NULL;
  double BestCost = NINFINITY;
  int Len;
  BackPoint *It = Path.Next;
  double LenPen,Add = 0.0;

  do {
    if (isinf(It->Additional)) break;
    Add += It->Additional;

    // when pen == 0 or StartPos is -1, this means we reach the
    // extremities of the sequence. Therefore we account for an
    // optimistic penality (MinPen) given that the length can only be
    // longer than the actual length. Because we use -1/Data_Len+1 as
    // indicators for sequence extremities we substract 1 from the
    // length in this case. Finally, in all cases, we discount the
    // Slope penalty on the length.

    if (pen && (It->StartPos >= 0)) {
      Len = (Forward ? pos-It->StartPos : It->StartPos-pos);
      LenPen = PenD[Len] - PenD.FinalSlope*Len;
    } else {
      Len = (Forward ? pos-It->StartPos-1 : It->StartPos-pos-1);
      LenPen = PenD.MinPen(Len) - PenD.FinalSlope*Len;
    }

    if ((Add + It->Cost - LenPen) > BestCost) {
      BestCost = Add+It->Cost - LenPen;
      BestBP = It;
    }
    if (It->Optimal && Len >= PenD.MaxLen) break;

    It = It -> Next;
  } while (It != Path.Next);
  

  *cost = BestCost;
  return BestBP;
}
// ----------------------------------------------------------------
// BackTrace and build a prediction object
// ----------------------------------------------------------------
Prediction* Track :: BackTrace (int MinCDSLen, int Forward)
{
  Prediction *pred = new Prediction();
  BackPoint  *It   = Path.Next;
  int  pos;
  char etat;
  int igpos = -1;
  int prevpos,prevstate;
  int CDSlen = 0;

  // initialisation by the terminal state
  pos  = (Forward ? It->StartPos : It->StartPos+1);
  etat = It->State;
  It   = It->Origin;
  pred->add(pos, etat);

  prevpos = pos;
  prevstate = etat;
  if (etat == InterGen) igpos = pos;

  while (It != NULL) {

    pos  = (Forward ? It->StartPos : It->StartPos+1);
    etat = It->State;
    It   = It->Origin;

    //    printf("pos %d etat %d CDSlen %d igpos %d prevpos %d\n",
    //	   pos,etat,CDSlen,igpos,prevpos);

    if (prevstate <= TermR3)  // codant
      CDSlen += prevpos-pos;

    prevpos = pos;
    prevstate = etat;

    if ((etat == InterGen) && (pos >=0))
      if (CDSlen <= MinCDSLen) {
	CDSlen = 0;
	pred->poptill(igpos);
	//	printf("CDSLen %d, je pope jusqu'a pos=%d (exclus)\n",
	//	       CDSlen,igpos);
      }
      else {
	igpos = pos;
	CDSlen = 0;
	pred->add(pos, etat);
      }
    else if (pos >=0)
      pred->add(pos, etat);
  }

  if (Forward) {
    pred->setPos(0, pred->getPos(0)-1);
  }
  else {
    pred->reversePred();
    pred->setPos(0,INITIALSHIFT+1-pred->getPos(0));
  }

  return pred;
}
// ----------------------------------------------------------------
// Dumps the contents of all the Backpoints in the path
// ----------------------------------------------------------------
void Track :: Dump ()
{
  BackPoint* It = Path.Next;

  printf("Number of BP allocated: %d\n\n",NumBP);

  do {
    printf("pos = %d, state = %d, cost = %f%s, additional = %f",
	   It->StartPos,It->State,It->Cost,
	   (It->Optimal ? "*" : ""), It->Additional);
    if (It->Origin)
      printf(" Or.state = %d Or.pos = %d",It->Origin->State,It->Origin->StartPos);
    printf("\n");
    It = It->Next;
  }  while (It != Path.Next);
}

