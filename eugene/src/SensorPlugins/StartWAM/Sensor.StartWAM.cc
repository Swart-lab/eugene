 //================================================================
//           Copyright (c) 2003 by INRA. All rights reserved.
//                 Redistribution is not permitted without
//                 the express written permission of INRA.
//                     Mail : tschiex@toulouse.inra.fr
//-----------------------------------------------------------------------------------------

// File                : EuGeneTk/SensorPlugins/StartWAM/Sensor.StartWAM.cc
// Description    : Definition of a start codon detection sensor based on a
//                         Weight Array Model
// Authors         : P.Bardou, S.Foissac, M.J.Cros, A.Moisan, T.Schiex

#include "Sensor.StartWAM.h"
#include <ctype.h>

extern Parameters PAR;

/*******************************************************
 **                  SensorStartWAM                   **
 *******************************************************/

// ----------------------
//  Default constructor.
// ----------------------
SensorStartWAM :: SensorStartWAM (int n) : Sensor(n)
{
  type = Type_Start;

  char modelfilename[FILENAME_MAX+1];
  MarkovianOrder= PAR.getI("StartWAM.MarkovianOrder");
  ScaleCoef = PAR.getD("StartWAM.ScaleCoef");
  ScalePenalty= PAR.getD("StartWAM.ScalePenalty");
  NbNtBeforeATG = PAR.getI("StartWAM.NbNtBeforeATG");
  NbNtAfterATG = PAR.getI("StartWAM.NbNtAfterATG");
  MotifLength= NbNtBeforeATG + 3 + NbNtAfterATG;
  strcpy(modelfilename,PAR.getC("EuGene.PluginsDir"));
  strcat(modelfilename,PAR.getC("StartWAM.modelfilename"));
  PlotScoreIncrease= 7.0;

  WAModel= new WAM(MarkovianOrder, MotifLength,"ACGT", modelfilename);
}

// ----------------------
//  Default destructor.
// ----------------------
SensorStartWAM :: ~SensorStartWAM ()
{
  delete WAModel;
}

// ---------------------
//  Init splice.
// ----------------------
void SensorStartWAM :: Init (DNASeq *X)
{
  if (PAR.getI("Output.graph")) Plot(X);
}

// ----------------
// Scaling
// ----------------
double SensorStartWAM :: ScaleWAMScore (double WAMScore) 
{
  return ( ScaleCoef  * WAMScore + ScalePenalty ) ;
}

// -----------------------
//  GiveInfo signal start.
// -----------------------
void SensorStartWAM :: GiveInfo (DNASeq *X, int pos, DATA *d)
{
  int i, j;
  double score;
  char* MotifExtnd = new char[MotifLength+MarkovianOrder+2]; 
  // Motif Extended = Motif + amount context 
  // (for markovian estimation of the first motif letter)
  MotifExtnd[MotifLength+MarkovianOrder+1] ='\0';

  ////////// START Forward (need enough context) //////////////
  if ( (X->IsEStart(pos,1)) &&
       (pos-NbNtBeforeATG-MarkovianOrder > 0) && 
       (pos+2+NbNtAfterATG < X->SeqLen) ) {
    score=0.0;
    j=0;
    for (i= pos-NbNtBeforeATG-MarkovianOrder; i<= pos+2+NbNtAfterATG; i++) {
      MotifExtnd[j]= toupper((*X)[i]);
      j++;
    }

    d->sig[DATA::Start].weight[Signal::Forward] += 
      ScaleWAMScore (WAModel->ScoreTheMotif(MotifExtnd));
  }

  ////////// START Reverse (need enough context)   //////////////

  if ( (X->IsEStart(pos-1,-1)) &&
       (pos-1+NbNtBeforeATG+MarkovianOrder < X->SeqLen) &&
       (pos-3-NbNtAfterATG > 0) ) {
    score=0.0;
    j=0;
    for (i= pos-1+NbNtBeforeATG+MarkovianOrder; i >= pos-3-NbNtAfterATG; i--) {
      MotifExtnd[j] = toupper((*X)(i));
      j++;
    }

    d->sig[DATA::Start].weight[Signal::Reverse] += 
      ScaleWAMScore (WAModel->ScoreTheMotif(MotifExtnd));
  }
}
// ----------------------------
//  Normalize Score for Plot
// ----------------------------
inline double SensorStartWAM :: NormalizePlot(double x, double n) 
{
  return Min(x,n)/n;
}

// ----------------------------
//  Plot Sensor information
// ----------------------------
void SensorStartWAM :: Plot(DNASeq *X)
{  
  int pos;
  double plotweight;
  DATA data;

  for(pos=0; pos < X-> SeqLen ; pos++){

    data.sig[DATA::Start].weight[Signal::Forward] = 0.0 ;
    data.sig[DATA::Start].weight[Signal::Reverse] = 0.0 ;

    GiveInfo (X, pos, &data);

    if (data.sig[DATA::Start].weight[Signal::Forward] != 0) {
      plotweight= data.sig[DATA::Start].weight[Signal::Forward] + PlotScoreIncrease;
      if (data.sig[DATA::Start].weight[Signal::Forward] > 0)
      	PlotBarF(pos,(pos%3)+1,0.5,NormalizePlot(plotweight,10.0),2);
    }
    if (data.sig[DATA::Start].weight[Signal::Reverse] != 0) {
      plotweight= data.sig[DATA::Start].weight[Signal::Reverse] + PlotScoreIncrease;
      if (data.sig[DATA::Start].weight[Signal::Reverse] > 0)
	PlotBarF(pos,-((X->SeqLen-pos)%3)-1,0.5,NormalizePlot(plotweight,10.0),2);
    }
  }
}

// ------------------
//  Post analyse
// ------------------
void SensorStartWAM :: PostAnalyse(Prediction *pred)
{
}
