//==============================================================
//           Copyright (c) 2003 by INRA. All rights reserved.
//                 Redistribution is not permitted without
//                 the express written permission of INRA.
//                     Mail : tschiex@toulouse.inra.fr
//-----------------------------------------------------------------------------------------

// File                : EuGeneTk/SensorPlugins/StartWAM/Sensor.StartWAM.h
// Description    : Definition of a start codon detection sensor based on a 
//                         Weight Array Model (WAM)
// Authors         : P.Bardou, S.Foissac, M.J.Cros, A.Moisan, T.Schiex

#ifndef  SENSOR_StartWAM_INCLUDED
#define  SENSOR_StartWAM_INCLUDED

#include "../../EuGene/Sensor.h"
#include "../0_SensorTk/WAM.h"

/*************************************************************
 **                 SensorStartWAM                    **
 *************************************************************/
class SensorStartWAM : public Sensor
{
 private:
  int NbNtBeforeATG;
  int NbNtAfterATG;
  int MotifLength;
  int MarkovianOrder; // order of the markov models in the WAM
  double ScaleCoef; // coefficient for the WAM score scaling
  double ScalePenalty; //  penality for the WAM score scaling
  double PlotScoreIncrease;
  WAM* WAModel;
  double ScaleWAMScore (double WAMScore);  //  (score= coef * WAMscore - pen)
  inline double NormalizePlot (double x, double n);  // normalization for the plot

 public:
  SensorStartWAM  (int n, DNASeq *X);
  virtual ~SensorStartWAM   ();
  virtual void Init       (DNASeq *);
  virtual void GiveInfo   (DNASeq *, int, DATA *);
  virtual void Plot       (DNASeq *);
  virtual void PostAnalyse(Prediction *);
};

extern "C" SensorStartWAM * builder0(int n, DNASeq *X) {  return new SensorStartWAM(n, X); }

#endif
