#include "Sensor.EuStop.h"

extern Parameters PAR;

/*************************************************************
 **                      SensorEuStop                       **
 *************************************************************/
// ----------------------
//  Default constructor.
// ----------------------
SensorEuStop :: SensorEuStop (int n) : Sensor(n)
{
}

// ----------------------
//  Default destructor.
// ----------------------
SensorEuStop :: ~SensorEuStop ()
{
}

// ----------------------
//  Init stop.
// ----------------------
void SensorEuStop :: Init (DNASeq *X)
{
  type = Type_Stop;
  stopP = PAR.getD("EuStop.stopP");
  
  if (PAR.getI("Output.graph")) Plot(X);
}

// -----------------------
//  GiveInfo signal stop.
// -----------------------
void SensorEuStop :: GiveInfo (DNASeq *X, int pos, DATA *d)
{
  double theStop;
  
  if ((theStop = X->IsStop(pos-3,1)) != 0.0) {
    d->sig[DATA::Stop].weight[Signal::Forward] += -stopP+log(theStop);
    d->sig[DATA::Stop].weight[Signal::ForwardNo] += log(1.0-theStop); 
  }
      
  if ((theStop = X->IsStop(pos+2,-1)) != 0.0) {
    d->sig[DATA::Stop].weight[Signal::Reverse] += -stopP+log(theStop);
    d->sig[DATA::Stop].weight[Signal::ReverseNo] += log(1.0-theStop); 
  }
}

// ----------------------------
//  Plot Sensor information
// ----------------------------
void SensorEuStop :: Plot(DNASeq *X)
{
  for(int i = 0;  i <= X->SeqLen;  i++) {
    if (X->IsStop(i-3,1)  != 0.0) 
      PlotBarF(i,(i%3)+1,0.1,0.2,1);

    if(X->IsStop(i+2,-1) != 0.0) 
      PlotBarF(i,-((X->SeqLen-i)%3)-1,0.1,0.2,1);
  }
}

// ------------------
//  Post analyse
// ------------------
void SensorEuStop :: PostAnalyse(Prediction *pred)
{
}
