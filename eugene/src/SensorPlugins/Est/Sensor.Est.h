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
// File:     Sensor.Est.h
// Contents: Sensor Est
// ------------------------------------------------------------------

#ifndef  SENSOR_EST_H_INCLUDED
#define  SENSOR_EST_H_INCLUDED

#include "../../Sensor.h"
#include "../../Hits.h"

/*************************************************************
 **                       SensorEst                         **
 *************************************************************/
class SensorEst : public Sensor
{
 private:
  std::vector<int>           vPos;
  std::vector<unsigned char> vESTMatch;
  std::vector<int>::iterator iter;
  int    index;
  unsigned char *ESTMatch;
  Hits   **HitTable;
  int    NumEST;
  double estP, utrP;
  double DonorThreshold;
  double spliceBoost;
  int    estM, utrM;
  int    ppNumber;
  int    stepid;
  int    N;
  
  Hits** ESTAnalyzer(FILE *, unsigned char *, int, int *, DNASeq *);
  void   ESTSupport (Prediction *pred,   int Tdebut,      int Tfin,
		     int debut, int fin, Hits **HitTable, int Size);
  void   FEASupport (Prediction *pred,   int Tdebut,      int Tfin,
		     int debut, int fin, Hits **HitTable, int Size, int NumG);
  int    LenSup     (Hits **HitTable, std::vector<int> vSupEstI,
		     int index, int beg, int end);

 public:
  SensorEst               (int n, DNASeq *X);
  virtual ~SensorEst      ();
  virtual void Init       (DNASeq *);
  virtual void GiveInfo   (DNASeq *, int, DATA *);
  virtual void Plot       (DNASeq *);
  virtual void PostAnalyse(Prediction *);
};

extern "C" SensorEst * builder0( int n, DNASeq *X) { return new SensorEst(n, X);}

#endif
