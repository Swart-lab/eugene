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
  
  indexR = indexF = 0;
  
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

  if (PAR.getI("Output.graph")) Plot(X);
}

// -----------------------
//  GiveInfo signal stop.
// -----------------------
void SensorEuStop :: GiveInfo (DNASeq *X, int pos, DATA *d)
{
  if((indexF != 0                 &&  vPosF[indexF-1] >= pos) ||
     (indexF < (int)vPosF.size()  &&  vPosF[indexF]   <  pos))
    {
      iter = lower_bound(vPosF.begin(), vPosF.end(), pos);
      if(*iter == pos) {
	d->Stop[0] += vValF[iter-vPosF.begin()];
	indexF = iter-vPosF.begin() + 1;
      }
      else indexF = iter-vPosF.begin();
    }
  else if( indexF < (int)vPosF.size()  &&  vPosF[indexF] == pos )
    {
      d->Stop[0] += vValF[indexF];
      indexF++;
    }

  if((indexR != 0                 &&  vPosR[indexR-1] >= pos) ||
     (indexR < (int)vPosR.size()  &&  vPosR[indexR]   <  pos))
    {
      iter = lower_bound(vPosR.begin(), vPosR.end(), pos);
      if(*iter == pos) {
	d->Stop[1] += vValR[iter-vPosR.begin()];
	indexR = iter-vPosR.begin() + 1;
      }
      else indexR = iter-vPosR.begin();
    }
  else if( indexR < (int)vPosR.size()  &&  vPosR[indexR] == pos )
    {
      d->Stop[1] += vValR[indexR];
      indexR++;
    }
}

// ----------------------------
//  Plot Sensor information
// ----------------------------
void SensorEuStop :: Plot(DNASeq *X)
{
  for (int i =0; i < (int)vPosF.size(); i++)
    PlotBarF(vPosF[i],(vPosF[i]%3)+1,0.1,0.2,1);

  for (int i =0; i < (int)vPosR.size(); i++)
    PlotBarF(vPosR[i],-((X->SeqLen-vPosR[i])%3)-1,0.1,0.2,1);
}

// ------------------
//  Post analyse
// ------------------
void SensorEuStop :: PostAnalyse(Prediction *pred)
{
}
