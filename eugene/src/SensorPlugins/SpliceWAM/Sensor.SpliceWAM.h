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
  double AccScaleCoef;
  double AccScalePenalty;
  double DonScaleCoef; // coef for the WAM score
  double DonScalePenalty; //  penalty for the WAM score (score= ScaleCoef * WAMscore - ScalePenalt
  static int NbNtBeforeGT, NbNtAfterGT, DonorSiteLength;
  static int NbNtBeforeAG, NbNtAfterAG, AcceptorSiteLength;
  static int MarkovianOrder;    // MarkovianOrder of the markov models in the Weight Array Model
  static WAM* DonWAModel;
  static WAM* AccWAModel;
  static bool IsInitialized;

 public:
  SensorSpliceWAM  (int n, DNASeq *X);
  virtual ~SensorSpliceWAM   ();
  virtual void Init       (DNASeq *);
  virtual void GiveInfo   (DNASeq *, int, DATA *);
  virtual void Plot       (DNASeq *);
  virtual void PostAnalyse(Prediction *);
};

extern "C" SensorSpliceWAM * builder0(int n, DNASeq *X) {  return new SensorSpliceWAM(n, X); }

int SensorSpliceWAM::NbNtBeforeGT, SensorSpliceWAM::NbNtAfterGT, SensorSpliceWAM::DonorSiteLength;
int SensorSpliceWAM::NbNtBeforeAG, SensorSpliceWAM::NbNtAfterAG, SensorSpliceWAM::AcceptorSiteLength;
int SensorSpliceWAM::MarkovianOrder;   
WAM* SensorSpliceWAM::DonWAModel;
WAM* SensorSpliceWAM::AccWAModel;
bool SensorSpliceWAM::IsInitialized = false;


#endif
