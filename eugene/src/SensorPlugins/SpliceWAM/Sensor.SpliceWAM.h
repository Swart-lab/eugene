//=================================================
//           Copyright (c) 2002 by INRA. All rights reserved.
//                 Redistribution is not permitted without
//                 the express written permission of INRA.
//                     Mail : tschiex@toulouse.inra.fr
//-----------------------------------------------------------------------------------------

// File                 : EuGeneTk/SensorPlugins/SpliceWAM/Sensor.SpliceWAM.h
// Description  : A simple splice site detection sensor based on a Weight Array Model
// Authors         : P.Bardou, S.Foissac, M.J.Cros, A.Moisan, T.Schiex

#ifndef  SENSOR_SPLICEWAM_INCLUDED
#define  SENSOR_SPLICEWAM_INCLUDED

#include "../../EuGene/Sensor.h"
#include "../0_SensorTk/WAM.h"

/*************************************************************
 **                 SensorSpliceWAM                    **
 *************************************************************/
class SensorSpliceWAM : public Sensor
{
 private:
  int NbNtBeforeGT, NbNtAfterGT, DonorSiteLength;
  int NbNtBeforeAG, NbNtAfterAG, AcceptorSiteLength;
  int MarkovianOrder;    // MarkovianOrder of the markov models in the Weight Array Model
  double DonScaleCoef; // coef for the WAM score
  double DonScalePenalty; //  penalty for the WAM score (score= ScaleCoef * WAMscore - ScalePenalt
  double AccScaleCoef;
  double AccScalePenalty;
  WAM* DonWAModel;
  WAM* AccWAModel;

 public:
  SensorSpliceWAM  (int);
  virtual ~SensorSpliceWAM   ();
  virtual void Init       (DNASeq *);
  virtual void GiveInfo   (DNASeq *, int, DATA *);
  virtual void Plot       (DNASeq *);
  virtual void PostAnalyse(Prediction *);
};

extern "C" SensorSpliceWAM * builder0(int n) {  return new SensorSpliceWAM(n); }

#endif
