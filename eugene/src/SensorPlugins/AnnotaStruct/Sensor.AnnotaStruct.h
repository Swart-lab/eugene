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

#ifndef  SENSOR_ANNOTASTRUCT_H_INCLUDED
#define  SENSOR_ANNOTASTRUCT_H_INCLUDED

#include <algorithm>
#include "../../Sensor.h"

/*************************************************************
 **                       Signals Object                    **
 *************************************************************/
class Signals
{
 private:
  
 public:
  int   pos;
  int   type;
  int   edge;
  char  score[20];
    
  Signals  ();
  Signals  (int, int, int, char*);
  ~Signals ();
  void PrintS ();
  bool operator < (int i) { if (pos < i) return true; else return false; }
};

/*************************************************************
 **                      Contents Object                    **
 *************************************************************/
class Contents
{
 private:
  
 public:
  int   start;
  int   end;
  int   type;
  float score;
    
  Contents  ();
  Contents  (int, int, int, float);
  ~Contents ();
  void PrintC ();
};

/*************************************************************
 **                      SensorAnnotaStruct
 *************************************************************/
class SensorAnnotaStruct : public Sensor
{
 private:
  std::vector <Signals*>  vSig;
  std::vector <Contents*> vCon;
  int   PosSigGiveInfo, PosConGiveInfo;
  int   iSig, iCon;
  char  *fileExt;
  char  startPAR[20], stopPAR[20],   accPAR[20];
  char  donPAR[20],   tStartPAR[20], tStopPAR[20];
  float exonPAR, intronPAR, cdsPAR;
  
  void ReadAnnotaStruct(char[FILENAME_MAX+1]);

 public:
  SensorAnnotaStruct          (int n, DNASeq *X);
  virtual ~SensorAnnotaStruct ();
  virtual void Init           (DNASeq *);
  virtual void GiveInfo       (DNASeq *X, int, DATA *);
  virtual void Plot           (DNASeq *X);
  virtual void PostAnalyse    (Prediction *);
};

extern "C" SensorAnnotaStruct* builder0( int n, DNASeq *X) { return new SensorAnnotaStruct(n, X);}

#endif
