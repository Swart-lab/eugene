// ------------------------------------------------------------------
// Copyright (C) 2004 INRA <eugene@ossau.toulouse.inra.fr>
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
// File:     Sensor.RibosomalFrameShift.cc
// Contents: Sensor RibosomalFrameShift (Constant Pattern)
// ------------------------------------------------------------------

#include <ctype.h>

#include "Sensor.RibosomalFrameShift.h"
#include "../../MSensor.h"

extern Parameters PAR;
extern MasterSensor *MS;

/*************************************************************
 **                       SensorRibosomalFrameShift
 *************************************************************/

// ----------------------
// Default constructor.
// ----------------------
SensorRibosomalFrameShift :: SensorRibosomalFrameShift (int n) : Sensor(n)
{
  
  type = Type_FS;
    
  newStatePos = PAR.getI("RibosomalFrameShift.newStatePos", GetNumber());
  patType     = PAR.getC("RibosomalFrameShift.type",        GetNumber());
  pattern     = PAR.getC("RibosomalFrameShift.pat",         GetNumber());
  requiredEstSupport = PAR.getI("RibosomalFrameShift.requiredEstSupport", GetNumber());
  patLen      = strlen(pattern);
  
  for(int i=0; i<patLen; i++)
    pattern[i] = tolower(pattern[i]);
  
  // Type initialisation
  if(!strcmp(patType, "insertion1")) 
  {
    sigTypeIndex = DATA::Ins1;      
  }
  else if(!strcmp(patType, "insertion2")) 
  {
    sigTypeIndex = DATA::Ins2;      
  }
  else if(!strcmp(patType, "insertion3")) 
  {
    sigTypeIndex = DATA::Ins3;      
  }
  else if(!strcmp(patType, "deletion1"))  
  {
    sigTypeIndex = DATA::Del1;
  }
    else if(!strcmp(patType, "deletion2"))  
  {
    sigTypeIndex = DATA::Del2;
  }
    else if(!strcmp(patType, "deletion3"))  
  {
    sigTypeIndex = DATA::Del3;
  }
  else 
  {
    fprintf(stderr, "Error \"%s\" undefined type in the parameter file"
	    " for RibosomalFrameShift.pat[%d]\n", patType, GetNumber());
    fprintf(stderr, "Type must be : insertion[123] or deletion[123].\n");
    exit(2);
  }
}

// ----------------------
//  Default destructor.
// ----------------------
SensorRibosomalFrameShift :: ~SensorRibosomalFrameShift ()
{
}
// ----------------------
//  Init.
// ----------------------
void SensorRibosomalFrameShift :: Init (DNASeq *X)
{
  patP   = PAR.getD("RibosomalFrameShift.patP*",  GetNumber());
  patPNo = PAR.getD("RibosomalFrameShift.patPNo*",GetNumber());

  if (PAR.getI("Output.graph")) Plot(X);
}

// ------------------------
//  GiveInfo.
// ------------------------
void SensorRibosomalFrameShift :: GiveInfo (DNASeq *X, int pos, DATA *d)
{
  int  i;
  bool isSigF = true;
  bool isSigR = true;

  // STRAND +
  for(i=0; i < patLen; i++)
  {
    if((*X)[pos+i-newStatePos+1] != pattern[i]) 
    {
      isSigF = false;
      break;
    }
  }
  // Now check there is a EST support
  if (isSigF && requiredEstSupport && !(d->EstMatch & HitForward) )
  {
    //cout << "Motif found but no EST support strand + position " << pos << "\n";
    isSigF = false;
  }
  
  
  // STRAND -
  for(i=patLen-1 ; i > -1; i--)
  {
    if((*X)(pos-i+newStatePos-2) != pattern[i]) 
    {
      isSigR = false;
      break;
    }
  }
  // Now check there is a EST support
  if ( isSigR && requiredEstSupport && !(d->EstMatch & HitReverse) )
  {
    //cout << "Motif found but no EST support strand - position " << pos << "\n";
    isSigR = false;
  }
  

  if (isSigF) 
  {
    //cout << "Motif found strand +  " << pos << "\n";
    d->sig[sigTypeIndex].weight[Signal::Forward]   += patP;
    d->sig[sigTypeIndex].weight[Signal::ForwardNo] += patPNo;
  }
  if (isSigR) 
  {
    //cout << "Motif found strand -  " << pos << "\n";
    d->sig[sigTypeIndex].weight[Signal::Reverse]   += patP;
    d->sig[sigTypeIndex].weight[Signal::ReverseNo] += patPNo;
  }
}

// ----------------------------
//  Plot Sensor information.
// ----------------------------
void SensorRibosomalFrameShift :: Plot(DNASeq *X)
{
}

// ------------------
//  Post analyse.
// ------------------
void SensorRibosomalFrameShift :: PostAnalyse(Prediction *pred, FILE *MINFO)
{
  // voir SensorFrameShift
}
