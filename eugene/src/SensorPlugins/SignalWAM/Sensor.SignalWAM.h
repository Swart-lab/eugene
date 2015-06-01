// ------------------------------------------------------------------
// Copyright (C) 2014 INRA <eugene-help@lists.mulcyber.toulouse.inra.fr>
//
// This program is open source; you can redistribute it and/or modify
// it under the terms of the Artistic License (see LICENSE file).
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
//
// You should have received a copy of Artistic License along with
// this program; if not, please see http://www.opensource.org
//
// $Id$
// ------------------------------------------------------------------
// File:     Sensor.SignalWAM.h
// Contents: Sensor SignalWAM
// A simple signal detection sensor based on a Weight Array Model
// ------------------------------------------------------------------


#ifndef  SENSOR_SIGNALWAM_H_INCLUDED
#define  SENSOR_SIGNALWAM_H_INCLUDED

#include "../../Sensor.h"
#include "../0_SensorTk/WAM.cc"

/*************************************************************
 **                 SensorSignalWAM                      **
 *************************************************************/
class SensorSignalWAM : public Sensor
{
 private:
  WAM* WAModel;
  int upstreamLen;   // nb of nt in the motif before the signal
  int downstreamLen; // nb of nt in the motif after the signal 
  int sigType;          // type of the signal DATA::Start, DATA::Acc or DATA::Don
  char*  pattern;       // pattern of the signal. Example : 'GT' for a donor signal
  
  double scaleCoef;     // coef for the WAM score
  double scalePenalty;  //  penalty for the WAM score (score= ScaleCoef * WAMscore - ScalePenalty
  
  std::vector<int>    vPosF, vPosR;
  std::vector<double> vScoreF, vScoreR;
  int indexF, indexR;
  int positionGiveInfo;
  
 int motifLength;       // length of the motif (signal+ flanking regions)
 int newStatePos;
 int scoringMode;       // if mode=0, score the full motif
                        // if mode=1, score the coding site of the motif only
 int scoringMotifStart; // start position in the motif to compute the score
 int scoringMotifEnd;   // end position in the motif to compute the score
 
 bool allowedNonCanSite; // true to search non canonic spliced site
                        
 void AssertPatternLength(char*, int, int, char*);
 void SearchSignal(DNASeq *X);
 
 public:
  SensorSignalWAM  (int n, DNASeq *X);
  virtual ~SensorSignalWAM  ();
  virtual void Init       (DNASeq *);
  virtual void GiveInfo   (DNASeq *, int, DATA *);
  virtual void Plot       (DNASeq *);
  virtual void PostAnalyse(Prediction *, FILE *);

};

extern "C" SensorSignalWAM * builder0(int n, DNASeq *X) {  return new SensorSignalWAM(n, X); }



#endif
