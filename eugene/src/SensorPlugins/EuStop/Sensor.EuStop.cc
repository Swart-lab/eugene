#include "Sensor.EuStop.h"

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
  vPosF.clear();
  vValF.clear();
  vPosR.clear();
  vValR.clear();
}

// ----------------------
//  Init stop.
// ----------------------
void SensorEuStop :: Init (DNASeq *X)
{
  type = Type_Stop;
  
  iterR = iterF = 0;
  
  vPosF.clear();
  vValF.clear();
  vPosR.clear();
  vValR.clear();
    
  for(int i = 0;  i <= X->SeqLen;  i++) {
    if(X->IsStop(i-3,1) != 0.0) {
      vPosF.push_back( i );
      vValF.push_back( X->IsStop(i-3,1) );
    }
    if(X->IsStop(i+2,-1) != 0.0) {
      vPosR.push_back( i );
      vValR.push_back( X->IsStop(i+2,-1) );
    }
  }
}

// -----------------------
//  ResetIter.
// -----------------------
void SensorEuStop :: ResetIter ()
{
  iterF = iterR = 0;
}

// -----------------------
//  GiveInfo signal stop.
// -----------------------
void SensorEuStop :: GiveInfo (DNASeq *X, int pos, DATA *d)
{
  if( iterF <= (int)vPosF.size()  &&  vPosF[iterF] == pos ) {
    d->Stop[0] += vValF[iterF];
    iterF++;
  }

  if( iterR <= (int)vPosR.size()  &&  vPosR[iterR] == pos ) {
    d->Stop[1] += vValR[iterR];
    iterR++;
  }
}

// -------------------------
//  GiveInfoAt signal stop.
// -------------------------
void SensorEuStop :: GiveInfoAt (DNASeq *X, int pos, DATA *d)
{
  iter = find(vPosF.begin(), vPosF.end(), pos);
  if(iter != vPosF.end())
    d->Stop[0] += vValF[iter-vPosF.begin()];
  
  iter = find(vPosR.begin(), vPosR.end(), pos);
  if(iter != vPosR.end())
    d->Stop[1] += vValR[iter-vPosR.begin()];
}
