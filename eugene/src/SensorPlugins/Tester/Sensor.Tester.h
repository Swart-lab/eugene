/*****************************************************************************/
/*             Copyright (c) 2002 by INRA. All rights reserved.              */
/*                 Redistribution is not permitted without                   */
/*                 the express written permission of INRA.                   */
/*                     Mail : tschiex@toulouse.inra.fr                       */
/*---------------------------------------------------------------------------*/
/* File         : EuGeneTk/SensorPlugins/Tester/Sensor.Tester.h              */
/* Description  : Sensor tester                                              */
/* Authors      : P.Bardou, S.Foissac, M.J.Cros, A.Moisan, T.Schiex          */
/* History      : May 2003         	   		                     */
/*****************************************************************************/

#ifndef  SENSOR_TESTER_H_INCLUDED
#define  SENSOR_TESTER_H_INCLUDED

#ifdef HAVE_CONFIG_H
#include "../../EuGene/config.h"
#endif

#ifdef HAVE_STRINGS_H
#include <strings.h>
#endif

#include "../../EuGene/Sensor.h"
#include "../../EuGene/Dll.h"
#include "../../EuGene/Prediction.h"

/*************************************************************
 **                     SensorTester                        **
 *************************************************************/
class SensorTester : public Sensor
{
 private:
  int  nbTest;
  char **source;
  char seqName[FILENAME_MAX+1];
  Sensor       **sensor;
  FILE         **fp;
  Prediction   * gene;
  
  void  ReadCoord(char[FILENAME_MAX+1]);
  char* SigType_TF(int, int, char **);
  char* State(int);

 public:
  SensorTester (int n, DNASeq *X);
  virtual ~SensorTester   ();
  virtual void Init       (DNASeq *);
  virtual void GiveInfo   (DNASeq *, int, DATA *);
  virtual void Plot       (DNASeq *);
  virtual void PostAnalyse(Prediction *);
};

extern "C" SensorTester * builder0(int n, DNASeq *X) { return new SensorTester(n, X); }

#endif
