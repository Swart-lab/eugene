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
  transStart = exp(-PAR.getD("Transcript.Start"));
  transStop = exp(-PAR.getD("Transcript.Stop"));

  if (PAR.getI("Output.graph")) Plot(X);
}

// -----------------------
//  GiveInfo signal start.
// -----------------------
void SensorTranscript :: GiveInfo (DNASeq *X, int pos, DATA *d)
{
  d->tStart[0] = transStart;
  d->tStart[1] = transStart;
  d->tStop[0] = transStop;
  d->tStop[1] = transStop;
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
