#include "Sensor.EuStop.h"

/*************************************************************
 **                      SensorEuStop                       **
 *************************************************************/
// ----------------------
//  Default constructor.
// ----------------------
SensorEuStop :: SensorEuStop (int n) : Sensor(n)
{
  *Stop = NULL;
}

// ----------------------
//  Default destructor.
// ----------------------
SensorEuStop :: ~SensorEuStop ()
{
  delete [] Stop[0];
  delete [] Stop[1];
}

// ----------------------
//  Init stop.
// ----------------------
void SensorEuStop :: Init (DNASeq *X)
{
  type = Type_Stop;

  if(*Stop != NULL) {
    delete [] Stop[0];
    delete [] Stop[1];
  }

  Stop[0] = new REAL[X->SeqLen+1];
  Stop[1] = new REAL[X->SeqLen+1];
  
  for  (int i = 0;  i <= X->SeqLen;  i++) {
    Stop[0][i] = X->IsStop(i-3,1);
    Stop[1][i] = X->IsStop(i+2,-1);
  }
}

// -----------------------
//  GiveInfo signal stop.
// -----------------------
void SensorEuStop :: GiveInfo (DNASeq *X, int pos, DATA *d)
{
  d->Stop[0] += Stop[0][pos];
  d->Stop[1] += Stop[1][pos];
}
