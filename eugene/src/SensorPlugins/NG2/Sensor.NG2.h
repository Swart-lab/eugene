/*****************************************************************************/
/*             Copyright (c) 2002 by INRA. All rights reserved.              */
/*                 Redistribution is not permitted without                   */
/*                 the express written permission of INRA.                   */
/*                     Mail : tschiex@toulouse.inra.fr                       */
/*---------------------------------------------------------------------------*/
/* File         : EuGeneTk/SensorPlugins/NG2/Sensor.NG2.h                    */
/* Description  : Sensor NetGene2                                            */
/* Authors      : P.Bardou, S.Foissac, M.J.Cros, A.Moisan, T.Schiex          */
/* History      : May 2003                                                   */
/*****************************************************************************/

#ifndef  SENSOR_NG2_H_INCLUDED
#define  SENSOR_NG2_H_INCLUDED

#include "../../EuGene/Sensor.h"

/*************************************************************
 **                       SensorNG2                         **
 *************************************************************/
class SensorNG2 : public Sensor
{
 private:
  int PositionGiveInfo;

  std::vector<int>    vPosAccF, vPosAccR, vPosDonF, vPosDonR;
  std::vector<double> vValAccF, vValAccR, vValDonF, vValDonR;

  int iAccF, iAccR, iDonF, iDonR;
  double accB, accP, donB, donP;
  
  void ReadNG2F(char[FILENAME_MAX+1], int);
  void ReadNG2R(char[FILENAME_MAX+1], int);

 public:
  SensorNG2  (int n, DNASeq *X);
  virtual ~SensorNG2      ();
  virtual void Init       (DNASeq *);
  virtual void GiveInfo   (DNASeq *X, int, DATA *);
  virtual void Plot       (DNASeq *X);
  virtual void PostAnalyse(Prediction *);
};

extern "C" SensorNG2 * builder0( int n, DNASeq *X) {  return new SensorNG2(n, X);}


#endif
