/*****************************************************************************/
/*             Copyright (c) 2002 by INRA. All rights reserved.              */
/*                 Redistribution is not permitted without                   */
/*                 the express written permission of INRA.                   */
/*                     Mail : tschiex@toulouse.inra.fr                       */
/*---------------------------------------------------------------------------*/
/* File         : EuGeneTk/SensorPlugins/Plotter/Sensor.Plotter.h            */
/* Description  : Sensor Plotter                                             */
/* Authors      : P.Bardou, S.Foissac, M.J.Cros, A.Moisan, T.Schiex          */
/* History      : May 2003                                                   */
/*****************************************************************************/

#ifndef  SENSOR_PLOTTER_H_INCLUDED
#define  SENSOR_PLOTTER_H_INCLUDED

#include "../../EuGene/Sensor.h"

/*************************************************************
 **                      SensorPlotter
 *************************************************************/
class SensorPlotter : public Sensor
{
 private:
  int window;
  int plotGC, plotGC3, plotAT;
  
 public:
  SensorPlotter          (int);
  virtual ~SensorPlotter ();
  virtual void Init       (DNASeq *);
  virtual void GiveInfo   (DNASeq *X, int, DATA *);
  virtual void Plot       (DNASeq *X);
  virtual void PostAnalyse(Prediction *);
};

extern "C" SensorPlotter* builder0( int n ) { return new SensorPlotter(n);}

#endif
