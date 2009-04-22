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
// File:     Sensor.NStretch.cc
// Contents: Sensor NStretch
// ------------------------------------------------------------------

#include "Sensor.NStretch.h"

extern Parameters PAR;

/*************************************************************
 **                      SensorNStretch                   **
 *************************************************************/

// ----------------------
//  Default constructor.
// ----------------------
SensorNStretch :: SensorNStretch (int n, DNASeq *X) : Sensor(n)
{
  type = Type_Content;
}

// ----------------------
//  Default destructor.
// ----------------------
SensorNStretch :: ~SensorNStretch ()
{

}

// ----------------------
//  Init Markov.insideNStretch
// ----------------------
void SensorNStretch :: Init (DNASeq *X)
{
	// récupérer le paramètre MaxLength et stretchPenalty
	//allouer un tableau insideNStretch de char
	// analyser la séquence pour remplir le char
}

// -----------------------
//  GiveInfo Content Markov.
// -----------------------
void SensorNStretch :: GiveInfo (DNASeq *X, int pos, DATA *d)
{
    if (insideNStretch[pos] != 0)
    {
      for(int i=0;i<6;i++)
        d->contents[i] -= stretchPenalty; //Exon

      d->contents[DATA::IntronF] -= stretchPenalty;
      d->contents[DATA::IntronR] -= stretchPenalty;
      d->contents[DATA::InterG] -= stretchPenalty;
      d->contents[DATA::UTR5F] -= stretchPenalty;
      d->contents[DATA::UTR5R] -= stretchPenalty;
      d->contents[DATA::UTR3F] -= stretchPenalty;
      d->contents[DATA::UTR3R] -= stretchPenalty;
      d->contents[DATA::IntronUTRF] -= stretchPenalty;
      d->contents[DATA::IntronUTRR] -= stretchPenalty;

    }
}
// ----------------------------
//  Plot Sensor information
// ----------------------------
void SensorNStretch :: Plot(DNASeq *TheSeq)
{
}

// ------------------
//  Post analyse
// ------------------
void SensorNStretch :: PostAnalyse(Prediction *pred, FILE *MINFO)
{
}
