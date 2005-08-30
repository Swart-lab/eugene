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
  StartPos = 0;
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
  if (cost > Optimal-PenD.GetDelta(abs(pos-OptPos)) && 
      cost > Path.Next->Cost+Path.Next->Additional-PenD.GetDelta(abs(pos-Path.Next->StartPos))) {
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
BackPoint *Track :: BestUsable(int pos, double *cost, int pen)
{
  BackPoint *BestBP = NULL;
  double BestCost = NINFINITY;
  int Len;
  BackPoint *It = Path.Next;
  double LenPen,Add = 0.0;

  do {
    if (isinf(It->Additional)) break;
    Add += It->Additional;
    Len = abs(pos-It->StartPos);

    // when pen == 0 or Origin is NULL, this means we reach the
    // extremities of the sequence. Therefore we account for an
    // optimistic penality (MinPen) given that the length can only be
    // longer than the actual length. In all cases, we discount the
    // Slope penalty on the length. We further discount one on extremities 
    // (because -1 and Data_Len+1 are used).

    if (pen && It->Origin) 
      LenPen = PenD[Len] - PenD.FinalSlope*Len;
    else 
      LenPen = PenD.MinPen(Len) - PenD.FinalSlope*(Len-1);

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
  std::vector <int>         vPos;
  std::vector <signed char> vState;
  Prediction *pred;
  BackPoint  *It;
  int  pos;
  char etat;
  int  ntopop = -1;
  int  prevpos, prevstate;
  int  CDSlen = 0;

  // put state back on correct transitions for backward predictions
  if (!Forward) {
    It = Path.Next;
    etat = It->State;
    It = It->Origin;
    
    while (It != NULL) {
      prevstate = etat;
      etat = It->State;
      It->State = prevstate;
      It = It->Origin;
    }
  }

  It = Path.Next;
  
  // initialisation by the terminal state
  pos  = It->StartPos;
  etat = It->State;
  It   = It->Origin;
  if (pos >= 0) {
    vPos.push_back  ( pos  );
    vState.push_back( etat );
  }

  prevpos = pos;
  prevstate = etat;
  if (etat == InterGen) ntopop=0;

  //  printf("pos %d etat %d CDSlen %d prevpos %d\n", pos,etat,CDSlen,prevpos);

  while (It != NULL) {

    pos  = It->StartPos;
    etat = It->State;
    It   = It->Origin;

    //printf("pos %d etat %d CDSlen %d prevpos %d\n", pos,etat,CDSlen,prevpos);

    // We count CDS. According to the direction the state to account
    // for is the current styate (forward) or the previous one
    // (backward).
    if ((Forward ? prevstate : etat) <= TermR3) {
      CDSlen += abs(prevpos-pos); // codant
      //      printf("Extra %d CDS, Len = %d\n",abs(prevpos-pos),CDSlen);
    }
    //  this is one more state change to count for possible removal
    if (ntopop >= 0) ntopop++;

    if ((etat == InterGen)) {
      if ((CDSlen > 0)  && (CDSlen <= MinCDSLen) && (ntopop >= 0)) {
	//	printf("CDSLen %d, je pope %d elements\n", CDSlen,ntopop-1);
	// we want to keep the same change state. In forward, we keep
	// the first IG state (pop n-1) and don't put back this one.
	if (Forward) {
	  for (int i = 0; i < ntopop-1; i++) {
	    vPos.pop_back();
	    vState.pop_back();
	  }
	}
	// in Backward, we pop n but insert this one.
	else {
	  for (int i = 0; i < ntopop; i++) {
	    vPos.pop_back();
	    vState.pop_back();
	  }
	  vPos.push_back  ( pos  );
	  vState.push_back( etat );
	}
      }
      else {
	vPos.push_back  ( pos  );
	vState.push_back( etat );
      }
      CDSlen = 0;
      ntopop = 0;
    }
    else if (pos >= 0) {
      vPos.push_back  ( pos  );
      vState.push_back( etat );
    }
    
    prevpos   = pos;
    prevstate = etat;
  }
  
  vPos[0] -=  1;  
  
  if (Forward) {
    reverse(vState.begin(), vState.end());
    reverse(vPos.begin(),   vPos.end());
  }
  
  pred = new Prediction(vPos, vState);
  vPos.clear();
  vState.clear();

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

