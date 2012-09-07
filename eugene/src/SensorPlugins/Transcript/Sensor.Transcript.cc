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
// File:     Sensor.Transcript.cc
// Contents: Sensor Transcript
// ------------------------------------------------------------------

#include "Sensor.Transcript.h"

extern Parameters PAR;

#define NORM(x,n) (((n)+(Max(-(n),x)))/(n))

/*************************************************************
 **                      SensorTranscript                   **
 *************************************************************/

// ----------------------
//  Default constructor.
// ----------------------
SensorTranscript :: SensorTranscript (int n, DNASeq *X) : Sensor(n)
{
  type = Type_Start;
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
  transStart     = PAR.getD("Transcript.Start*",         GetNumber());
  transStop      = PAR.getD("Transcript.Stop*",          GetNumber());
  transStartNpc  = PAR.getD("Transcript.StartNpc*",      GetNumber());
  transStopNpc   = PAR.getD("Transcript.StopNpc*",       GetNumber());
  affectedStrand = PAR.getI("Transcript.AffectedStrand", GetNumber());

  if (PAR.getI("Output.graph")) Plot(X);
}

// -----------------------
//  GiveInfo signal start.
// -----------------------
void SensorTranscript :: GiveInfo (DNASeq *X, int pos, DATA *d)
{
  // Strand +
  if (affectedStrand != -1)
  {
	d->sig[DATA::tStart].weight[Signal::Forward]    -= transStart;
    d->sig[DATA::tStop].weight[Signal::Forward]     -= transStop;	
    d->sig[DATA::tStartNpc].weight[Signal::Forward] -= transStartNpc;
    d->sig[DATA::tStopNpc].weight[Signal::Forward]  -= transStopNpc;
  }
  // Strand -
  if (affectedStrand != 1)
  {
	d->sig[DATA::tStart].weight[Signal::Reverse]    -= transStart;
    d->sig[DATA::tStop].weight[Signal::Reverse]     -= transStop;	
    d->sig[DATA::tStartNpc].weight[Signal::Reverse] -= transStartNpc;
    d->sig[DATA::tStopNpc].weight[Signal::Reverse]  -= transStopNpc;
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
void SensorTranscript :: PostAnalyse(Prediction *pred, FILE *MINFO)
{
}
