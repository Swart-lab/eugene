#include "Sensor.Transcript.h"
#define NORM(x,n) (((n)+(Max(-(n),x)))/(n))

/*************************************************************
 **                      SensorTranscript                   **
 *************************************************************/

extern Parameters PAR;

// ----------------------
//  Default constructor.
// ----------------------
SensorTranscript :: SensorTranscript (int n) : Sensor(n)
{
}

// ----------------------
//  Default destructor.
// ----------------------
SensorTranscript :: ~SensorTranscript ()
{
}

// ----------------------
//  Init start.
// ----------------------
void SensorTranscript :: Init (DNASeq *X)
{
  type = Type_Start;
  transStart = PAR.getD("Transcript.Start");
  transStop = PAR.getD("Transcript.Stop");

  if (PAR.getI("Output.graph")) Plot(X);
}

// -----------------------
//  GiveInfo signal start.
// -----------------------
void SensorTranscript :: GiveInfo (DNASeq *X, int pos, DATA *d)
{
  for (int i =0; i <= Signal::Reverse; i++) {
    d->sig[DATA::tStart].weight[i] -= transStart;
    d->sig[DATA::tStop].weight[i] -= transStop;
  }
}

// ----------------------------
//  Plot Sensor information
// ----------------------------
void SensorTranscript :: Plot(DNASeq *X)
{
}

// ------------------
//  Post analyse
// ------------------
void SensorTranscript :: PostAnalyse(Prediction *pred)
{
}