/*****************************************************************************/
/*             Copyright (c) 2004 by INRA. All rights reserved.              */
/*                 Redistribution is not permitted  without                  */
/*                 the express written permission of  INRA.                  */
/*                   Mail : eugene@ossau.toulouse.inra.fr                    */
/*---------------------------------------------------------------------------*/
/* File         : EuGeneTk/SensorPlugins/Est/Sensor.Est.h                    */
/* Description  : Sensor Est                                                 */
/* Authors      : P.Bardou, S.Foissac, M.J.Cros, A.Moisan, T.Schiex          */
/* History      : July 2004                                                  */
/*****************************************************************************/

#ifndef  SENSOR_EST_H_INCLUDED
#define  SENSOR_EST_H_INCLUDED

#include "../../EuGene/Sensor.h"
#include "../../EuGene/Hits.h"

/*************************************************************
 **                       SensorEst                         **
 *************************************************************/
class SensorEst : public Sensor
{
 private:
  std::vector<int>           vPos;
  std::vector<unsigned char> vESTMatch;
  std::vector<int>::iterator iter;
  int    index;
  unsigned char *ESTMatch;
  Hits   **HitTable;
  int    NumEST;
  double estP, utrP;
  double DonorThreshold;
  double spliceBoost;
  int    estM, utrM;
  int    ppNumber;
  int    stepid;
  int    N;
  
  Hits** ESTAnalyzer(FILE *, unsigned char *, int, int *, DNASeq *);
  void   ESTSupport (Prediction *pred,   int Tdebut,      int Tfin,
		     int debut, int fin, Hits **HitTable, int Size);
  void   FEASupport (Prediction *pred,   int Tdebut,      int Tfin,
		     int debut, int fin, Hits **HitTable, int Size, int NumG);
  int    LenSup     (Hits **HitTable, std::vector<int> vSupEstI,
		     int index, int beg, int end);

 public:
  SensorEst               (int n, DNASeq *X);
  virtual ~SensorEst      ();
  virtual void Init       (DNASeq *);
  virtual void GiveInfo   (DNASeq *, int, DATA *);
  virtual void Plot       (DNASeq *);
  virtual void PostAnalyse(Prediction *);
};

extern "C" SensorEst * builder0( int n, DNASeq *X) { return new SensorEst(n, X);}

#endif
