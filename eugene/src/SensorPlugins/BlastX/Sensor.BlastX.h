/*****************************************************************************/
/*             Copyright (c) 2004 by INRA. All rights reserved.              */
/*                 Redistribution is not permitted  without                  */
/*                 the express written permission of  INRA.                  */
/*                   Mail : eugene@ossau.toulouse.inra.fr                    */
/*---------------------------------------------------------------------------*/
/* File         : EuGeneTk/SensorPlugins/BlastX/Sensor.BlastX.h              */
/* Description  : Sensor BlastX                                              */
/* Authors      : P.Bardou, S.Foissac, M.J.Cros, A.Moisan, T.Schiex          */
/* History      : July 2004                                                  */
/*****************************************************************************/

#ifndef  SENSOR_BLASTX_H_INCLUDED
#define  SENSOR_BLASTX_H_INCLUDED

#include "../../EuGene/Sensor.h"
#include "../../EuGene/Hits.h"

/*************************************************************
 **                     SensorBlastX                        **
 *************************************************************/
class SensorBlastX : public Sensor
{
 private:
  double *ProtMatch, *ProtMatchLevel;
  int    *ProtMatchPhase;
  std::vector<int>    vPos,     vPMPhase;
  std::vector<double> vPMLevel, vPMatch;
  std::vector<int>::iterator iter;
  char *levels;
  int    index;
  double keyBXLevel[10];
  int    minIn;
  int    blastxM;
  int    ppNumber;
  int    N;

  void LoadContentScore (DNASeq *);
  char ph06             (char);

  // For postprocess 2
  Hits **HitTable;
  int  NumProt;
  void ProtSupport (Prediction *pred, int debut, int fin,
		    Hits **HitTable,  int Size,  int NumG);
  int  LenSup      (Hits **HitTable, Prediction *pred,
		    std::vector<int> vSupProtI,
		    int index, int beg, int end);
  
 public:
  SensorBlastX  (int n, DNASeq *X);
  virtual ~SensorBlastX   ();
  virtual void Init       (DNASeq *);
  virtual void GiveInfo   (DNASeq *, int, DATA *);
  virtual void Plot       (DNASeq *);
  virtual void PostAnalyse(Prediction *);
};

extern "C" SensorBlastX * builder0( int n, DNASeq *X) { return new SensorBlastX(n, X); }
#endif
