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
// File:     Sensor.Template.h
// Contents: Sensor Template
// ------------------------------------------------------------------

#ifndef  SENSOR_TEMPLATE_INCLUDED
#define  SENSOR_TEMPLATE_INCLUDED

#include "../../Sensor.h"

/*************************************************************
 **                     SensorTemplate                    **
 *************************************************************/
class SensorTemplate : public Sensor
{
 private:

 public:
  SensorTemplate  (int n, DNASeq *X);
  virtual ~SensorTemplate   ();
  virtual void Init       (DNASeq *);
  virtual void GiveInfo   (DNASeq *, int, DATA *);
  virtual void Plot       (DNASeq *);
  virtual void PostAnalyse(Prediction *, FILE *);
};

extern "C" SensorTemplate * builder0(int n, DNASeq *X) {  return new SensorTemplate(n, X); }

#endif
