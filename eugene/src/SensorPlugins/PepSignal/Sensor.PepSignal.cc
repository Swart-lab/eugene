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
// File:     Sensor.PepSignal.cc
// Contents: Sensor Peptide Signal
// ------------------------------------------------------------------

#include "Sensor.PepSignal.h"
#define NORM(x,n) (((n)+(Max(-(n),x)))/(n))

/*************************************************************
 **                       SensorPepSignal                   **
 *************************************************************/

extern Parameters PAR;

// ----------------------
// Default constructor.
// ----------------------
SensorPepSignal :: SensorPepSignal (int n, DNASeq *X) : Sensor(n)
{
char *seqname;
char tempname[FILENAME_MAX+1];

	seqname = PAR.getC("fstname");

	fprintf(stderr, "Probing PepSignal (starts)....");  
	fflush(stderr);

	strcpy(tempname,seqname);
	strcat(tempname,".psignal");

	ReadPepSignalStarts(tempname, X->SeqLen);
	fprintf(stderr,"done\n");

	CheckStart(X,vPosF, vPosR);
}


// --------------------------
//  Read start file.
// --------------------------
void SensorPepSignal :: ReadPepSignalStarts(char *name, int Len)
{
  
  FILE *fp;
  char type[10];  // max len for "start_rev"
  int len = 0;
  int i,pos,end =0;
  double force;

  if (!(fp = fopen(name, "r"))) {
    fprintf(stderr, "Cannot open start file %s\n", name);
    exit(2);
  }

  len =1;
  while (!end) {
    i = fscanf(fp,"%d %s %lf %*s\n", &pos, type, &force);
    len ++;
    if (i == EOF) {end = 1; break;}
    if (i < 3) {end =2; break;}

    if (strcmp(type,"start") == 0) {
      vPosF.push_back(pos-1);
      vValF.push_back(force);
    } else if (strcmp(type,"start_rev") == 0) {
      vPosR.push_back(pos);
      vValR.push_back(force);
    } else end = 2;
  }
  fclose(fp);
  
  if (end ==2) {
    fprintf(stderr, "Error in  PepSignal start file %s, line %d\n", name, len);
    exit(2);
  }
}




// ----------------------
//  Default destructor.
// ----------------------
SensorPepSignal :: ~SensorPepSignal ()
{
  vPosF.clear();
  vValF.clear();
  vPosR.clear();
  vValR.clear();
}

// ----------------------
//  Init start.
// ----------------------
void SensorPepSignal :: Init (DNASeq *X)
{
  startP = PAR.getD("PepSignal.startP*",GetNumber());
  startB = PAR.getD("PepSignal.startB*",GetNumber());

  indexR = indexF = 0;
  PositionGiveInfo = -1;
  
  if (PAR.getI("Output.graph")) Plot(X);
}

// ------------------------
//  GiveInfo signal start.
// ------------------------
void SensorPepSignal :: GiveInfo (DNASeq *X, int pos, DATA *d)
{
  bool update = false;
  double f;

  if ( (PositionGiveInfo == -1) || (pos != PositionGiveInfo+1) ) update = true; // update indexes on vectors
  PositionGiveInfo = pos;
  
  // Start Forward
  if(!vPosF.empty()) {
    if (update) 
      indexF = lower_bound(vPosF.begin(), vPosF.end(), pos)-vPosF.begin();
    
    if((indexF<(int)vPosF.size()) && (vPosF[indexF] == pos)) {
      f = pow(vValF[indexF], startB)*(exp(-startP));
      d->sig[DATA::Start].weight[Signal::Forward] += log(f);
      d->sig[DATA::Start].weight[Signal::ForwardNo] += log(1.0-f);
      indexF++;
    }
  }
  
  // Start Reverse
  if (!vPosR.empty()) {
    if (update) 
      indexR = lower_bound(vPosR.begin(), vPosR.end(), pos)-vPosR.begin();
    
    if((indexR<(int)vPosR.size()) && (vPosR[indexR] == pos)) {
      f = pow(vValR[indexR], startB)*(exp(-startP));
      d->sig[DATA::Start].weight[Signal::Reverse] += log(f);
      d->sig[DATA::Start].weight[Signal::ReverseNo] += log(1.0-f);
      indexR++;
    }
  }


}

// ----------------------------
//  Plot Sensor information
// ----------------------------
void SensorPepSignal :: Plot(DNASeq *X)
{
  double f;

  for (int i =0; i < (int)vPosF.size(); i++) {
    f = pow(vValF[i], startB)*(exp(-startP));
    PlotBarF(vPosF[i], (vPosF[i]%3)+1, 0.5, NORM(log(f),4.0), 2);
  }

  for (int i =0; i < (int)vPosR.size(); i++) {
    f = pow(vValR[i], startB)*(exp(-startP));
    PlotBarF(vPosR[i], -((X->SeqLen-vPosR[i])%3)-1, 0.5, NORM(log(f),4.0), 2);
  }
}

// ------------------
//  Post analyse
// ------------------
void SensorPepSignal :: PostAnalyse(Prediction *pred)
{
}