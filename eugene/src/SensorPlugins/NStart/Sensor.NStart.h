/*****************************************************************************/
/*             Copyright (c) 2002 by INRA. All rights reserved.              */
/*                 Redistribution is not permitted without                   */
/*                 the express written permission of INRA.                   */
/*                     Mail : tschiex@toulouse.inra.fr                       */
/*---------------------------------------------------------------------------*/
/* File         : EuGeneTk/SensorPlugins/NStart/Sensor.NStart.cc             */
/* Description  : Sensor Netstart                                            */
/* Authors      : P.Bardou, S.Foissac, M.J.Cros, A.Moisan, T.Schiex          */
/* History      : May 2003                                                   */
/*****************************************************************************/

#ifndef  SENSOR_NSTART_H_INCLUDED
#define  SENSOR_NSTART_H_INCLUDED

#include "../../EuGene/Sensor.h"

/*************************************************************
 **                    SensorNStart                         **
 *************************************************************/
class SensorNStart : public Sensor
{
 private:
  int PositionGiveInfo;

  std::vector<int>    vPosF, vPosR;
  std::vector<double> vValF, vValR;

  int indexF, indexR;
  double startP, startB;
  
  void ReadNStartF (char *, int);
  void ReadNStartR (char *, int);

 public:
  SensorNStart   (int n, DNASeq *X);
  virtual ~SensorNStart   ();
  virtual void Init       (DNASeq *);
  virtual void GiveInfo   (DNASeq *X, int, DATA *);
  virtual void Plot       (DNASeq *X);
  virtual void PostAnalyse(Prediction *);
};

extern "C" SensorNStart* builder0( int n, DNASeq *X) { return new SensorNStart(n, X);}

#endif
