/*****************************************************************************/
/*             Copyright (c) 2004 by INRA. All rights reserved.              */
/*                 Redistribution is not permitted without                   */
/*                 the express written permission of INRA.                   */
/*                     Mail : tschiex@toulouse.inra.fr                       */
/*---------------------------------------------------------------------------*/
/* File         : EuGeneTk/SensorPlugins/PatConst/Sensor.PatConst.h          */
/* Description  : Sensor pattern const                                       */
/* Authors      : P.Bardou, S.Foissac, M.J.Cros, A.Moisan, T.Schiex          */
/* History      : Fev 2004         	   		                     */
/*****************************************************************************/

#ifndef  SENSOR_PATCONST_H_INCLUDED
#define  SENSOR_PATCONST_H_INCLUDED

#include "../../EuGene/Sensor.h"

/*************************************************************
 **                      SensorPatConst
 *************************************************************/
class SensorPatConst : public Sensor
{
 private:
  double patP;
  double patPNo;
  char*  pattern;
  char*  patType;
  int    newStatePos;
  int    sigTypeIndex;
  int    patLen;
  
 public:
  SensorPatConst          (int);
  virtual ~SensorPatConst ();
  virtual void Init       (DNASeq *);
  virtual void GiveInfo   (DNASeq *X, int, DATA *);
  virtual void Plot       (DNASeq *X);
  virtual void PostAnalyse(Prediction *);
};

extern "C" SensorPatConst* builder0( int n ) { return new SensorPatConst(n);}

#endif
