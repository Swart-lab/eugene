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
// File:     Sensor.h
// Contents: class Sensor
// ------------------------------------------------------------------

#ifndef  SENSOR_H_INCLUDED
#define  SENSOR_H_INCLUDED

#include <vector>
#include <algorithm>
#include <string.h>
#ifdef HAVE_STRINGS_H
#include <strings.h>
#endif

#include "Param.h"
#include "DNASeq.h"
#include "Prediction.h"
#include "Const.h"
#include "SensorIF.h"
extern "C" {
#include "GDIF/gdIF.h"
}

/*************************************************************
 **                        Sensor                           **
 *************************************************************/
class Sensor
{
 private:
  int instanceNumber;

 protected:
  void CheckStart   (DNASeq *,std::vector<int> a, std::vector<int> b);
  void CheckStop    (DNASeq *,std::vector<int> a, std::vector<int> b);
  void CheckSplices (DNASeq *,std::vector<int> a, std::vector<int> b, std::vector<int> c, std::vector<int> d);
  int  GetNumber    (void) { return instanceNumber; }

 public:
  unsigned char type;
  
  Sensor  (int n);
  virtual ~Sensor (void);
  virtual void Init       (DNASeq *X) = 0;
  virtual void GiveInfo   (DNASeq *, int, DATA *) = 0;
  virtual void Plot       (DNASeq *) = 0;
  virtual void PostAnalyse(Prediction *) = 0;
};

#endif

