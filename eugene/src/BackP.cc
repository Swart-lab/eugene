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
			    double cost, BackPoint *Or)
{
  BackPoint *It = this->Next;

  if (cost > NINFINITY) {
    It =  new BackPoint(state,pos,cost);
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
// BackTrace along a BackPoint and store the parsing in Choix
// ----------------------------------------------------------------
void BackPoint :: BackTrace (char *Choix)
{
  BackPoint* It = this->Next;
  int pos,i;
  char etat;
  
  pos = It->StartPos;
  etat =It->State;
  It  = It->Origin;
  
  do {
    for (i = pos; i > It->StartPos; i--) Choix[i] = etat;
    pos = It->StartPos;
    etat = It->State;
    It = It->Origin;
  }  while (It != NULL);
}

// ----------------------------------------------------------------
// Updates the cost 
// ----------------------------------------------------------------
void BackPoint :: Update(double cost)
{
  Next->Additional += cost;
}

// ----------------------------------------------------------------
// Returns the best BackPoint. If the SwitchType matches the mask
// given, then the BackPoint MUST BE at least len nuc. far from pos
// ----------------------------------------------------------------
BackPoint *BackPoint :: BestUsable(int pos, unsigned char mask, int len, REAL *cost)
{
  BackPoint *It = this->Next;
  double add = 0.0;

  do {
    add += It->Additional;
    if ((It->StartPos <= pos-len) || !(It->SwitchType & mask) || (It->StartPos < 0)) {
      *cost = add+It->Cost;
      return It;
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
    if ((It->StartPos <= pos-len) || (It->StartPos < 0)) {
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
