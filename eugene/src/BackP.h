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
// File:     BackP.h
// Contents: Definitions for a linear time shortest path with constraint alg.
// Comments: Une liste doublement chainee et circulaire de point de retour
// arriere pour le plus court chemin avec contrainte de duree minimum.
// Le premier element Path (dans Track) ne doit pas etre
// efface. Path.Next est le dernier insere, Path.Prev est le plus vieux
// ??? A t'on vraiment besoin de Prev ?
// ------------------------------------------------------------------

#ifndef  BACKP_H_INCLUDED
#define  BACKP_H_INCLUDED


#include "Const.h"
#include "Prediction.h"
#include "PenaltyDist.h"

class BackPoint
{
  friend class Track;
  // private:
 public:
  signed char State;
  char Optimal;
  int StartPos;  
  double Cost;
  double Additional;
  BackPoint *Origin;
  
 public:
  
  BackPoint *Next;
  BackPoint *Prev;    
  
  BackPoint ();
  BackPoint  (char state, int pos, double cost);
  ~BackPoint();

  inline void InitState(int state, int start) { State = state; StartPos = start;};
  void Print();
};


class Track 
{
 public: 

  int NumBP;
  BackPoint Path;
  PenaltyDist PenD;
  double Optimal;
  int OptPos;

  inline void LoadPenalty(char* name) { PenD.LoadPenaltyDist(name); };
  inline void Update(double cost) {  Path.Next->Additional += cost; Optimal += cost;};
  inline void PayTheSlope() {Path.Next->Additional -= PenD.FinalSlope;
                             Optimal  -= PenD.FinalSlope;  };
  void InsertNew(char state,  int pos, double cost, BackPoint *Or);
  void ForceNew(char state, int pos, double cost, BackPoint *Or);
  BackPoint *BestUsable(int pos, double *cost, int pen = 1);
  Prediction* BackTrace(int MinCDSLen, int Forward = 1);
  void Dump();
  void Zap();

  Track ();
  ~Track();
  
};

#endif
