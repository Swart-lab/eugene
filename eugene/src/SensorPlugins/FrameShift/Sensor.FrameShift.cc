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
// File:     Sensor.FrameShift.cc
// Contents: Sensor FrameShift
// ------------------------------------------------------------------

#include "Sensor.FrameShift.h"

extern Parameters PAR;

#define NORM(x,n) (((n)+(Max(-(n),x)))/(n))

/*************************************************************
 **                      SensorFrameShift                   **
 *************************************************************/

// ----------------------
//  Default constructor.
// ----------------------
SensorFrameShift :: SensorFrameShift (int n, DNASeq *X) : Sensor(n)
{
    type = Type_FS;
}

// ----------------------
//  Default destructor.
// ----------------------
SensorFrameShift :: ~SensorFrameShift ()
{
}

// ----------------------
//  Init start.
// ----------------------
void SensorFrameShift :: Init (DNASeq *X)
{
    insProb = -(PAR.getD("FrameShift.Ins*"));
    delProb = -(PAR.getD("FrameShift.Del*"));

    if (PAR.getI("Output.graph")) Plot(X);
}

// -----------------------
//  GiveInfo frameshift
// -----------------------
void SensorFrameShift :: GiveInfo (DNASeq *X, int pos, DATA *d)
{   
    d->sig[DATA::Ins1].weight[Signal::Forward] = insProb;
    d->sig[DATA::Ins2].weight[Signal::Forward] = insProb;
    d->sig[DATA::Ins3].weight[Signal::Forward] = insProb;
    d->sig[DATA::Ins1].weight[Signal::Reverse] = insProb;
    d->sig[DATA::Ins2].weight[Signal::Reverse] = insProb;
    d->sig[DATA::Ins3].weight[Signal::Reverse] = insProb;
    
    d->sig[DATA::Del1].weight[Signal::Forward] = delProb;
    d->sig[DATA::Del2].weight[Signal::Forward] = delProb;
    d->sig[DATA::Del3].weight[Signal::Forward] = delProb;
    d->sig[DATA::Del1].weight[Signal::Reverse] = delProb;
    d->sig[DATA::Del2].weight[Signal::Reverse] = delProb;
    d->sig[DATA::Del3].weight[Signal::Reverse] = delProb;
}

// ----------------------------
//  Plot Sensor information
// ----------------------------
void SensorFrameShift :: Plot(DNASeq *X)
{
}

// ------------------
//  Post analyse
// ------------------
void SensorFrameShift :: PostAnalyse(Prediction *pred, FILE *MINFO)
{
    short int frame;
    short int frameBack = 0;
    int posFs = -1;

    if (PAR.getI("Output.graph"))
    {
        for (int i=0; i<pred->nbGene; i++)
        {
            for (int j=0; j<pred->vGene[i]->nbFea(); j++)
            {
                if (pred->vGene[i]->vFea[j]->IsCodingExon())
                {
                    frame = pred->vGene[i]->vFea[j]->frame;
                    if (posFs != -1)  // Frameshift plot
                    {
                        PlotLine(posFs, posFs, frame, frameBack, 0.4, 0.4, 1);
                    }
                    posFs     = pred->vGene[i]->vFea[j]->end;
                    frameBack = frame;
                }
                else                // Ig, UTR, or Intron
                    posFs = -1;
            }
        }
    }
}
