/*****************************************************************************/
/*             Copyright (c) 2002 by INRA. All rights reserved.              */
/*                 Redistribution is not permitted without                   */
/*                 the express written permission of INRA.                   */
/*                     Mail : tschiex@toulouse.inra.fr                       */
/*---------------------------------------------------------------------------*/
/* File         : EuGeneTk/SensorPlugins/SPred/Sensor.SPred.h                */
/* Description  : Sensor SplicerPredictor                                    */
/* Authors      : P.Bardou, S.Foissac, M.J.Cros, A.Moisan, T.Schiex          */
/* History      : May 2003                                                   */
/*****************************************************************************/

#ifndef  SENSOR_SPred_H_INCLUDED
#define  SENSOR_SPred_H_INCLUDED

#include "../../EuGene/Sensor.h"

/*************************************************************
 **                     SensorSPred                         **
 *************************************************************/
class SensorSPred : public Sensor
{
 private:
  int PositionGiveInfo;

  std::vector<int>    vPosAccF, vPosAccR, vPosDonF, vPosDonR;
  std::vector<double> vValAccF, vValAccR, vValDonF, vValDonR;

  int iAccF, iAccR, iDonF, iDonR;
  double accP, accB, donP, donB;
  
  void ReadSPredF(char[FILENAME_MAX+1], int);
  void ReadSPredR(char[FILENAME_MAX+1], int);
  
 public:
  SensorSPred   (int n, DNASeq *X);
  virtual ~SensorSPred    ();
  virtual void Init       (DNASeq *);
  virtual void GiveInfo   (DNASeq *X, int, DATA *);
  virtual void Plot       (DNASeq *X);
  virtual void PostAnalyse(Prediction *);
};

extern "C" SensorSPred * builder0( int n, DNASeq *X) {  return new SensorSPred(n, X);}

#endif
