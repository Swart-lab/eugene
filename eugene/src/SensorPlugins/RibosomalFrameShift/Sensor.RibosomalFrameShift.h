// ------------------------------------------------------------------
// Copyright (C) 2012 INRA <eugene@ossau.toulouse.inra.fr>
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

// ------------------------------------------------------------------
// File:     Sensor.RibosomalFrameShift.h
// Contents: Sensor RibosomalFrameShift (Constant Pattern)
// ------------------------------------------------------------------

#ifndef  SENSOR_RIBOSOMALFRAMESHIFT_H_INCLUDED
#define  SENSOR_RIBOSOMALFRAMESHIFT_H_INCLUDED

#include "../../Sensor.h"

/*************************************************************
 **                      SensorRibosomalFrameShift
 *************************************************************/
class SensorRibosomalFrameShift : public Sensor
{
 private:
  double patP;
  double patPNo;
  char*  pattern;
  char*  patType;
  int    newStatePos;
  int    sigTypeIndex;
  int    patLen;
  int    requiredEstSupport;
  
 public:
  SensorRibosomalFrameShift          (int);
  virtual ~SensorRibosomalFrameShift ();
  virtual void Init       (DNASeq *);
  virtual void GiveInfo   (DNASeq *X, int, DATA *);
  virtual void Plot       (DNASeq *X);
  virtual void PostAnalyse(Prediction *, FILE *);
};

extern "C" SensorRibosomalFrameShift* builder0( int n ) { return new SensorRibosomalFrameShift(n);}

#endif
