/*****************************************************************************/
/*             Copyright (c) 2002 by INRA. All rights reserved.              */
/*                 Redistribution is not permitted without                   */
/*                 the express written permission of INRA.                   */
/*                     Mail : tschiex@toulouse.inra.fr                       */
/*---------------------------------------------------------------------------*/
/* File         : EuGeneTk/SensorPlugins/ATGpr/Sensor.ATGpr.h                */
/* Description  : Sensor ATGpr                                               */
/* Authors      : P.Bardou, S.Foissac, M.J.Cros, A.Moisan, T.Schiex          */
/* History      : May 2003                                                   */
/*****************************************************************************/

#ifndef  SENSOR_ATGPR_H_INCLUDED
#define  SENSOR_ATGPR_H_INCLUDED

#include "../../EuGene/Sensor.h"

/*************************************************************
 **                     SensorATGpr                         **
 *************************************************************/
class SensorATGpr : public Sensor
{
 private:
  int PositionGiveInfo;
  std::vector<int>  vPosF, vPosR;
  std::vector<REAL> vValF, vValR;
  int indexF, indexR;
  double startP, startB;
  
  void ReadATGprF (char[FILENAME_MAX+1], int);
  void ReadATGprR (char[FILENAME_MAX+1], int);

 public:
  SensorATGpr   (int);
  virtual ~SensorATGpr   ();
  virtual void Init       (DNASeq *);
  virtual void GiveInfo   (DNASeq *X, int, DATA *);
  virtual void Plot       (DNASeq *X);
  virtual void PostAnalyse(Prediction *);
};

extern "C" SensorATGpr* builder0( int n ) { return new SensorATGpr(n);}

#endif
