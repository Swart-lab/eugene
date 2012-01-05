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
// File:     Sensor.MarkovConst.cc
// Contents: Sensor MarkovConst
// ------------------------------------------------------------------

#include "Sensor.MarkovConst.h"

extern Parameters PAR;

/*************************************************************
 **                      SensorMarkovConst                   **
 *************************************************************/

// ----------------------
//  Default constructor.
// ----------------------
SensorMarkovConst :: SensorMarkovConst (int n, DNASeq *X) : Sensor(n)
{
  type = Type_Content;
}

// ----------------------
//  Default destructor.
// ----------------------
SensorMarkovConst :: ~SensorMarkovConst ()
{
}

// ----------------------
//  Init Markov.
// ----------------------
void SensorMarkovConst :: Init (DNASeq *X)
{
  transCodant  = PAR.getD("MarkovConst.Coding*", GetNumber()); //Exon
  transIntron  = PAR.getD("MarkovConst.Intron*", GetNumber()); //Intron
  transIntronU = PAR.getD("MarkovConst.IntronUTR*", GetNumber()); //IntronUTR
  transInter   = PAR.getD("MarkovConst.Inter*", GetNumber());  //InterG
  transUTR5    = PAR.getD("MarkovConst.UTR5*", GetNumber());   //UTR5
  transUTR3    = PAR.getD("MarkovConst.UTR3*", GetNumber());   //UTR3
  transUIR     = PAR.getD("MarkovConst.UIR*", GetNumber());   //UIR
  transRNA     = PAR.getD("MarkovConst.RNA*", GetNumber());   // RNA

  minGC = PAR.getD("MarkovConst.minGC",GetNumber());
  maxGC = PAR.getD("MarkovConst.maxGC",GetNumber());

  affectedStrand = PAR.getI("MarkovConst.affectedStrand",GetNumber());
}

// -----------------------
//  GiveInfo Content Markov.
// -----------------------
void SensorMarkovConst :: GiveInfo (DNASeq *X, int pos, DATA *d)
{
    if ((X->Markov0[BitG] + X->Markov0[BitC]) > minGC &&
            (X->Markov0[BitG] + X->Markov0[BitC]) <= maxGC)
    {
        if (affectedStrand != -1)
        {
            for (int i=0; i<3; i++)
                d->contents[i] += log(transCodant); //Exon
            d->contents[DATA::IntronF]    += log(transIntron);
            d->contents[DATA::UTR5F]      += log(transUTR5);
            d->contents[DATA::UTR3F]      += log(transUTR3);
            d->contents[DATA::IntronUTRF] += log(transIntronU);
            d->contents[DATA::RNAF]       += log(transRNA);
            d->contents[DATA::UIRF]       += log(transUIR);
        }
        if (affectedStrand != 1)
        {
            for (int i=3; i<6 ; i++)
                d->contents[i] += log(transCodant); //Exon
            d->contents[DATA::IntronR]    += log(transIntron);
            d->contents[DATA::UTR5R]      += log(transUTR5);
            d->contents[DATA::UTR3R]      += log(transUTR3);
            d->contents[DATA::IntronUTRR] += log(transIntronU);
            d->contents[DATA::RNAR]       += log(transRNA);
            d->contents[DATA::UIRR]       += log(transUIR);
        }
        d->contents[DATA::InterG] += log(transInter);
    }
}

// ----------------------------
//  Plot Sensor information
// ----------------------------
void SensorMarkovConst :: Plot(DNASeq *TheSeq)
{
}

// ------------------
//  Post analyse
// ------------------
void SensorMarkovConst :: PostAnalyse(Prediction *pred, FILE *MINFO)
{
}
