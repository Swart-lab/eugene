/*****************************************************************************/
/*             Copyright (c) 2002 by INRA. All rights reserved.              */
/*                 Redistribution is not permitted without                   */
/*                 the express written permission of INRA.                   */
/*                     Mail : tschiex@toulouse.inra.fr                       */
/*---------------------------------------------------------------------------*/
/* File         : EuGeneTk/SensorPlugins/SMachine/Sensor.SMachine.h                    */
/* Description  : Sensor NetGene2                                            */
/* Authors      : P.Bardou, S.Foissac, M.J.Cros, A.Moisan, T.Schiex          */
/* History      : Mars 2004                                                   */
/*****************************************************************************/

#ifndef  SENSOR_SMACHINE_H_INCLUDED
#define  SENSOR_SMACHINE_H_INCLUDED

#include "../../EuGene/Sensor.h"


/*************************************************************
 **                       SensorSMachine                         **
 *************************************************************/
class SensorSMachine : public Sensor
{
 private:
  int PositionGiveInfo;

  std::vector<int>    vPosAccF, vPosAccR, vPosDonF, vPosDonR;
  std::vector<double> vValAccF, vValAccR, vValDonF, vValDonR;

  std::vector<int>    vPosF, vPosR;
  std::vector<double> vValF, vValR;

  int iAccF, iAccR, iDonF, iDonR;
  double accB, accP, donB, donP;

  int indexF, indexR;
  double startP, startB;

  void ReadSMachineSplices(char *, int);
  void ReadSMachineStarts(char *, int);

  void SpliceMachine();
  
 public:
  SensorSMachine  (int n, DNASeq *X);
  virtual ~SensorSMachine ();
  virtual void Init       (DNASeq *);
  virtual void GiveInfo   (DNASeq *X, int, DATA *);
  virtual void Plot       (DNASeq *X);
  virtual void PostAnalyse(Prediction *);
};

extern "C" SensorSMachine * builder0( int n, DNASeq *X) {  return new SensorSMachine(n, X);}


#endif
