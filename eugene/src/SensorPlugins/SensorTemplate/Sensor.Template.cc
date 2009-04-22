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
// File:     Sensor.Template.cc
// Contents: Sensor Template
// ------------------------------------------------------------------

#include "Sensor.Template.h"

extern Parameters PAR;

/*************************************************************
 **                      SensorTemplate                   **
 *************************************************************/

// ----------------------
//  Default constructor.
// ----------------------
SensorTemplate :: SensorTemplate (int n, DNASeq *X) : Sensor(n)
{
  type = Type_Content;
}

// ----------------------
//  Default destructor.
// ----------------------
SensorTemplate :: ~SensorTemplate ()
{
}

// ----------------------
//  Init Markov.
// ----------------------
void SensorTemplate :: Init (DNASeq *X)
{

}

// -----------------------
//  GiveInfo Content Markov.
// -----------------------
void SensorTemplate :: GiveInfo (DNASeq *X, int pos, DATA *d)
{

}

// ----------------------------
//  Plot Sensor information
// ----------------------------
void SensorTemplate :: Plot(DNASeq *TheSeq)
{
}

// ------------------
//  Post analyse
// ------------------
void SensorTemplate :: PostAnalyse(Prediction *pred, FILE *MINFO)
{
}
