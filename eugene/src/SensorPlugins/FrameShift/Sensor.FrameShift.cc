#include "Sensor.FrameShift.h"
#define NORM(x,n) (((n)+(Max(-(n),x)))/(n))

/*************************************************************
 **                      SensorFrameShift                   **
 *************************************************************/

extern Parameters PAR;

// ----------------------
//  Default constructor.
// ----------------------
SensorFrameShift :: SensorFrameShift (int n) : Sensor(n)
{
}

// ----------------------
//  Default destructor.
// ----------------------
SensorFrameShift :: ~SensorFrameShift ()
{
}

// ----------------------
//  Init start.
// ----------------------
void SensorFrameShift :: Init (DNASeq *X)
{
  type = Type_Start;
  insProb = -exp(PAR.getD("FrameShift.Ins"));
  delProb = -exp(PAR.getD("FrameShift.Del"));
  if (PAR.getI("Output.graph")) Plot(X);
}

// -----------------------
//  GiveInfo frameshift
// -----------------------
void SensorFrameShift :: GiveInfo (DNASeq *X, int pos, DATA *d)
{
    d->sig[DATA::Ins].weight[Signal::Forward] = insProb;
    d->sig[DATA::Ins].weight[Signal::Reverse] = insProb;
    d->sig[DATA::Del].weight[Signal::Forward] = delProb;
    d->sig[DATA::Del].weight[Signal::Reverse] = delProb;
}

// ----------------------------
//  Plot Sensor information
// ----------------------------
void SensorFrameShift :: Plot(DNASeq *X)
{
}

// ------------------
//  Post analyse
// ------------------
void SensorFrameShift :: PostAnalyse(Prediction *pred)
{
  int state, stateBack;
  int posFs = -1;
  
  if (PAR.getI("Output.graph")) {
    for(int i=pred->size()-1; i!=-1; i--) {
      state = pred->getState(i);
      if(state <= ExonR3)
	if(posFs != -1) { // Frameshift plot
	  if(state <= ExonF3)	 
	    PlotLine(posFs, posFs, state+1, stateBack+1, 0.4, 0.4, 1);
	  else
	    PlotLine(posFs, posFs, -state+2, -stateBack+2, 0.4, 0.4, 1);
	  
	  posFs = pred->getPos(i);
	  stateBack = state;
	}
	else {            // No frameshift
	  posFs = pred->getPos(i);
	  stateBack = state;
	}
      else                // Ig, UTR, or Intron 
	posFs = -1;
    }
  }
}
