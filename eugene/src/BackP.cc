//
//   T. Schiex
//
//     File:  BackP.cc
//     Version:  1.0
//
//    Copyright (c) 2000 by Thomas Schiex All rights reserved.
//    Redistribution is not permitted without the express written
//    permission of the authors.
//
//  Definitions for a linear time shortest path with constraint alg.
//

#include "BackP.h"
#include <stdio.h>

// ----------------------------------------------------------------
// Default constructor.
// ----------------------------------------------------------------
BackPoint :: BackPoint  ()
  {
    Optimal = true;
    State = -1;
    StartPos = -1;
    SwitchType = SwitchAny;
    Cost = Additional = 0.0;
    Next = Prev = Origin = NULL;
  }

// ----------------------------------------------------------------
//  Default constructor.
// ----------------------------------------------------------------
BackPoint :: BackPoint  (char state, int pos, double cost)
{
  Optimal = true;
  State = state;
  StartPos = pos;
  SwitchType = SwitchAny;
  Cost = cost;
  Additional = 0.0;
  Next = Prev = Origin = NULL;
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
// Insert  a new backpoint
// ----------------------------------------------------------------
void BackPoint :: InsertNew(char state, unsigned char Switch, int pos, 
			    double cost, BackPoint *Or, bool opt)
{
  BackPoint *It = this->Next;

  if (cost > NINFINITY) {
    It =  new BackPoint(state,pos,cost);
    It->Optimal = opt;
    It->Next = this->Next;
    It->Prev = this->Next->Prev;
    It->Origin = Or;
    It->SwitchType = Switch;
    this->Next->Prev = It;
    this->Next = It;
  }
}

// ----------------------------------------------------------------
// Prints the BackPoint contents
// ----------------------------------------------------------------
void BackPoint :: Print()
{
  printf("pos = %d, state = %d\n", StartPos,State);
}

// ----------------------------------------------------------------
// Dumps the contents of all the BackPoints
// ----------------------------------------------------------------
void BackPoint :: Dump ()
{
  BackPoint* It = this->Next;
  double add = 0.0;
    
  do {
    add += It->Additional;
    printf("pos = %d, state = %d, cost = %f, real cost = %f\n",
	   It->StartPos,It->State,It->Cost,add+It->Cost);
    It = It->Next;
  }  while (It != this->Next);
}

// ----------------------------------------------------------------
// BackTrace along a BackPoint and build a prediction object
// ----------------------------------------------------------------
Prediction* BackPoint :: BackTrace ()
{
  Prediction *pred = new Prediction();
  BackPoint  *It   = this->Next;
  int  pos;
  char etat;
  
  pos  = It->StartPos;
  etat = It->State;
  It   = It->Origin;
  
  while (It != NULL) {
    pred->add(pos, etat);
    pos  = It->StartPos;
    etat = It->State;
    It   = It->Origin;
  }
  pred->setPos(0, pred->getPos(0)-1);

  return pred;
}
// ----------------------------------------------------------------
// Returns the best BackPoint. If the SwitchType matches the mask
// given, then the BackPoint MUST BE at least len nuc. far from pos
// else it must be at least len2 nuc. far
// ----------------------------------------------------------------
BackPoint *BackPoint :: BestUsable(int pos, unsigned char mask, int len, REAL *cost, int len2)
{
  BackPoint *It = this->Next;
  BackPoint *BestBP = NULL;
  double BestScore = NINFINITY;
  double add = 0.0;

  do {
    add += It->Additional;
    // la condition filtre les aiguillages interdits: trop pres (len2) ou trop pres (len) et de mauvais type
    if ((It->StartPos <= pos-len) || ( (!(It->SwitchType & mask))&& (It->StartPos <= pos-len2)) || (It->StartPos < 0)) {
      // on regarde si on ameliore
      if ((add+It->Cost) > BestScore) {
	BestScore = add+It->Cost;
	BestBP = It;
      }
      // si c'est un aiguillage autorise optimal, c'est fini sinon on continue.
      if (It->Optimal) {
	*cost = BestScore;
	return BestBP;
      }
    }
    It = It -> Next;
  } while (It != this->Next);

  *cost = NINFINITY;
  return NULL;
}
// ----------------------------------------------------------------
// Returns the best BackPoint and the BackPoint is at least len
// nuc. far from pos
// ----------------------------------------------------------------
BackPoint *BackPoint :: StrictBestUsable(int pos, int len, REAL *cost)
{
  BackPoint *It = this->Next;
  double add = 0.0;

  do {
    add += It->Additional;
    if (It->Optimal && ((It->StartPos <= pos-len) || (It->StartPos < 0))) {
      *cost = add+It->Cost;
      return It;
    }
    It = It -> Next;
  } while (It != this->Next);

  *cost = NINFINITY;
  return NULL;
}

// ----------------------------------------------------------------
// Zap  a whole list
// ----------------------------------------------------------------
void BackPoint :: Zap()
{
  while (Next != this) delete this->Next;
  delete this;
}
