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
// File:     Sensor.Repeat.cc
// Contents: Sensor Repeat
// ------------------------------------------------------------------

#include "Sensor.Repeat.h"

extern Parameters PAR;

/*************************************************************
 **                        SensorRepeat                     **
 *************************************************************/

// ----------------------
// Default constructor.
// ----------------------
SensorRepeat :: SensorRepeat (int n, DNASeq *X) : Sensor(n)
{
  char tempname[FILENAME_MAX+1];
  FILE* ncfile;
  int deb, fin;

  type = Type_Content;

  fprintf(stderr,"Reading Intergenic regions....");
  fflush(stderr);
  
  strcpy(tempname,PAR.getC("fstname"));
  strcat(tempname,".ig");
  ncfile = FileOpen(NULL,tempname, "r");
  
  while (fscanf(ncfile,"%d %d\n", &deb, &fin) != EOF)  {
    int s = (int)vDeb.size();
    deb   = Max(1,deb)-1;
    fin   = Min(X->SeqLen,fin)-1;
    if(deb > fin  ||
       (s != 0  &&  vFin[s-1] >= deb)) {
      fprintf(stderr,"\nError in ig file %s, line %d\n", tempname, s+1);
      exit(2);
    }
    vDeb.push_back( deb );
    vFin.push_back( fin );
  }
  fprintf(stderr,"done\n");
}

// ----------------------
//  Default destructor.
// ----------------------
SensorRepeat :: ~SensorRepeat ()
{
  vDeb.clear();
  vFin.clear();
}

// ----------------------
//  Init Repeat.
// ----------------------
void SensorRepeat :: Init (DNASeq *X)
{
  UTRPenalty = PAR.getD("Repeat.UTRPenalty*");
  intronPenalty = PAR.getD("Repeat.IntronPenalty*");
  exonPenalty = PAR.getD("Repeat.ExonPenalty*");
  
  index = 0;
  PositionGiveInfo = -1;

  if (PAR.getI("Output.graph")) Plot(X);
}

// --------------------------
//  GiveInfo Content Repeat.
// --------------------------
void SensorRepeat :: GiveInfo (DNASeq *X, int pos, DATA *d)
{
  bool update = false;

  if ( (PositionGiveInfo == -1) || (pos != PositionGiveInfo+1) ) update = true; // update indexes on vectors
  PositionGiveInfo = pos;

  if(!vFin.empty()) {
    if (update) 
      index = lower_bound(vFin.begin(), vFin.end(), pos) - vFin.begin();

    if (index < (int) vFin.size())
      if (pos >= vDeb[index]) { 
	for(int i=DATA::ExonF1; i<=DATA::ExonR3; i++)	 // Exon(6)
	  d->contents[i] -= exonPenalty;
	for(int i=DATA::IntronF; i<=DATA::IntronR; i++)  // Intron(2)
	  d->contents[i] -= intronPenalty; 
	for(int i=DATA::UTR5F; i<=DATA::IntronUTRR; i++) // UTR(4)+IntronUTR(2)
	  d->contents[i] -= UTRPenalty;
 	if (pos == vFin[index]) index++;
      }
  }
}

// ----------------------------
//  Plot Sensor information
// ----------------------------
void SensorRepeat :: Plot(DNASeq *TheSeq)
{
  int i;
  
  for (i =0; i < (int)vDeb.size(); i++)
    PlotRepeat(vDeb[i],vFin[i]);
}

// ------------------
//  Post analyse
// ------------------
void SensorRepeat :: PostAnalyse(Prediction *pred)
{
}
