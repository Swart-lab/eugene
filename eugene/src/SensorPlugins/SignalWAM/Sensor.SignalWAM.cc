// ------------------------------------------------------------------
// Copyright (C) 2014 INRA <eugene-help@lists.mulcyber.toulouse.inra.fr>
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
// File:     Sensor.SignalWAM.cc
// Contents: Sensor SignalWAM
// A signal detection sensor based on a Weight Array Model
// ------------------------------------------------------------------

#include <ctype.h>
#include "Sensor.SignalWAM.h"

extern Parameters PAR;

#define NORM(x,n) ( (Min(x,n))/(n))
#define plotscoreincrease 10

/*************************************************************
 **                  SensorSignalWAM                   **
 *************************************************************/

// ----------------------
//  Default constructor.
// ----------------------
SensorSignalWAM :: SensorSignalWAM(int n, DNASeq *X) : Sensor(n)
{
    char* signalType;
    signalType       = PAR.getC("SignalWAM.type",          GetNumber()); // start, acc or don
    pattern          = PAR.getC("SignalWAM.pat",           GetNumber()); // 'ATG' for instance
    markovianOrder   = PAR.getI("SignalWAM.markovOrder",   GetNumber()); // 0 1 2
    upstreamLen      = PAR.getI("SignalWAM.uc",            GetNumber());
    downstreamLen    = PAR.getI("SignalWAM.dc",            GetNumber());
    scoringMode      = PAR.getI("SignalWAM.scoringMode",   GetNumber());
    
    
    char modelfilename[FILENAME_MAX + 1];
    strcpy(modelfilename, PAR.getC("SignalWAM.filePrefix", GetNumber()));
    sigLength    = strlen(pattern);
    motifLength  = upstreamLen + sigLength + downstreamLen;
    
    scoringMotifStart = 1;
    scoringMotifEnd   = motifLength;

    if (!strcmp(signalType, "start"))
    {
        AssertPatternLength(pattern, 3, GetNumber(), "start");
        sigType     = DATA::Start;
        type        = Type_Start;
        newStatePos = 1;

        if (scoringMode == 1) //codingside
        {
            scoringMotifStart = upstreamLen+4;
            scoringMotifEnd   = motifLength;
        }
    }
    else if (!strcmp(signalType, "donor"))
    {
        AssertPatternLength(pattern, 2, GetNumber(), "donor");
        sigType     = DATA::Don;
        type        = Type_Don;
        newStatePos = 1;
        
        if (scoringMode == 1) //codingside
        {
            scoringMotifStart = 1;
            scoringMotifEnd   = upstreamLen;
        }
    }
    else if (!strcmp(signalType, "acceptor"))
    {
        AssertPatternLength(pattern, 2, GetNumber(), "acceptor");
        sigType     = DATA::Acc;
        type        = Type_Acc;
        newStatePos = 3;
        if (scoringMode == 1) //codingside
        {
            scoringMotifStart = upstreamLen+3;
            scoringMotifEnd   = motifLength;
        }
    }
    else
    {
        fprintf(stderr, "Error \"%s\" undefined type in the parameter file"
                        " for SignalWAM.type[%d]\n", signalType, GetNumber());
        fprintf(stderr, "Type must be : start, acceptor or donor.\n");
        exit(2);  
    }

    WAModel = new WAM(markovianOrder, motifLength, "ACGT", modelfilename);
}

// ----------------------
//  Default destructor.
// ----------------------
SensorSignalWAM :: ~SensorSignalWAM()
{
    delete [] WAModel;
}

// ---------------------
//  Init
// ----------------------
void SensorSignalWAM :: Init(DNASeq *X)
{
    scaleCoef    = PAR.getD("SignalWAM.scaleCoef*", GetNumber());
    scalePenalty = PAR.getD("SignalWAM.scalePenalty*", GetNumber());

    if (PAR.getI("Output.graph"))
        Plot(X);
}

// ---------------------
//  AssertPatternLength
// ----------------------
void SensorSignalWAM :: AssertPatternLength(char* pattern, int expectedlength, int instNb, char* sigType)
{
    int sigLength = strlen(pattern);

    if (sigLength != expectedlength)
    {
        fprintf(stderr, "Error: %s is not allowed for SignalWAM.pat[%d] (expected %s signal length is %d)\n", pattern, instNb,  sigType, expectedlength);
        exit(2);
    }
}

// -----------------------
//  GiveInfo.
// -----------------------
void SensorSignalWAM :: GiveInfo(DNASeq *X, int pos, DATA *d)
{
    int i, j;
    double score;
    char* Site = new char[motifLength + 1];
    Site[motifLength] = '\0';

    int sigStartPos;
    int sigStopPos;

    // -----------------------------------------------------------------
    // Forward strand
    // -----------------------------------------------------------------
    bool isSigF = true;
    sigStartPos = pos - newStatePos + 1;
    sigStopPos  = sigStartPos + sigLength - 1;

    // check there is a signal at the current position on forward strand
    for (i = 0; i < sigLength; i++)
    {
        if ((*X)[sigStartPos + i] != tolower(pattern[i]))
        {
            isSigF = false;
            break;
        }
    }

    // If signal detects, extract the context sequence and compute the score
    if ((isSigF) && (sigStartPos - upstreamLen > 0) && (sigStopPos + downstreamLen < X->SeqLen))
    {
        score = 0.0;
        j     = 0;

        for (i = sigStartPos - upstreamLen; i <= sigStopPos + downstreamLen; i++)
        {
            Site[j] = toupper((*X)[i]);
            j++;
        }
        //if (scoringMode == 2)
        //       score = WAModel->ScoreTheSubMotif(Site, 1, upstreamLen) * WAModel->ScoreTheSubMotif(Site, upstreamLen+sigLength+1, motifLength);

        score = WAModel->ScoreTheMotif(Site, scoringMotifStart, scoringMotifEnd);
       
        //fprintf(stdout, "SITE POSITION = %d : mode=%d, [%d-%d], site=%s, score = %f\n", pos, scoringMode, scoringMotifStart, scoringMotifEnd, Site, score);
        d->sig[sigType].weight[Signal::Forward] += (scaleCoef * score) + scalePenalty;
    }

    // -----------------------------------------------------------------
    // Reverse strand
    // ------------------------------------------------------------------
    bool isSigR = true;
    sigStartPos = pos + newStatePos - 2;
    sigStopPos  = sigStartPos - sigLength + 1;

    // check there is a signal at the current position on forward strand
    for (i = sigLength - 1; i > -1; i--)
    {
        if ((*X)(sigStartPos - i) != tolower(pattern[i]))
        {
            isSigR = false;
            break;
        }
    }

    // If signal detects, extract the context sequence and compute the score
    if ((isSigR) && (sigStartPos + upstreamLen < X->SeqLen) && (sigStopPos - downstreamLen > 0))
    {
        score = 0.0;
        j     = 0;
        
        for (i = sigStartPos + upstreamLen; i >= sigStopPos - downstreamLen; i--)
        {
            Site[j] = toupper((*X)(i));
            j++;
        }

        //if (scoringMode == 2)
        //score = WAModel->ScoreTheSubMotif(Site, 1, upstreamLen) * WAModel->ScoreTheSubMotif(Site, upstreamLen+sigLength+1, motifLength);
        score = WAModel->ScoreTheMotif(Site, scoringMotifStart, scoringMotifEnd);
                
        //fprintf(stdout, "REV SITE POSITION = %d : mode=%i, [%d-%d], site=%s, score = %f\n", pos, scoringMode, scoringMotifStart, scoringMotifEnd, Site, score);
        d->sig[sigType].weight[Signal::Reverse] += (scaleCoef * score) + scalePenalty;
    }

    delete [] Site;
}


// ----------------------------
//  Plot Sensor information
// ----------------------------
void SensorSignalWAM :: Plot(DNASeq *X)
{
    int i, pos;
    DATA data;

    for (pos = 0; pos < X-> SeqLen ; pos++)
    {

        data.sig[sigType].weight[Signal::Forward] = data.sig[sigType].weight[Signal::Reverse] = 0.0;

        GiveInfo(X, pos, &data);

        if (data.sig[sigType].weight[Signal::Forward] != 0)
        {
            data.sig[sigType].weight[Signal::Forward] += plotscoreincrease;

            if (data.sig[sigType].weight[Signal::Forward] > 0)
            {
                if (sigType == DATA::Acc)
                    PlotAcc(pos, 1, NORM(data.sig[sigType].weight[Signal::Forward], 20.0));
                else
                    if (sigType == DATA::Don)
                        PlotDon(pos, 1, NORM(data.sig[sigType].weight[Signal::Forward], 20.0));
                    else
                        if (sigType == DATA::Start)
                            PlotStart(pos, (pos % 3) + 1, NORM(data.sig[sigType].weight[Signal::Forward], 10.0));
            }
        }

        if (data.sig[sigType].weight[Signal::Reverse] != 0)
        {
            data.sig[sigType].weight[Signal::Reverse] += plotscoreincrease;

            if (data.sig[sigType].weight[Signal::Reverse] > 0)
            {
                if (sigType == DATA::Acc)
                    PlotAcc(pos, -1, NORM(data.sig[sigType].weight[Signal::Reverse], 20.0));
                else
                    if (sigType == DATA::Don)
                        PlotDon(pos, -1, NORM(data.sig[sigType].weight[Signal::Reverse], 20.0));
                    else
                        if (sigType == DATA::Start)
                            PlotStart(pos, -((X->SeqLen - pos) % 3) - 1, NORM(data.sig[sigType].weight[Signal::Reverse], 10.0));
            }
        }

        /* to print all the sites: (differents colors if negatives)
            if (data.sig[DATA::Acc].weight[Signal::Forward] != 0 )
              PlotBarF(pos,4,0.5,NORM((fabs)(data.sig[DATA::Acc].weight[Signal::Forward]),20.0), (data.sig[DATA::Acc].weight[Signal::Forward] >0) ? 4 : 7 );
            if (data.sig[DATA::Don].weight[Signal::Forward] != 0 )
              PlotBarF(pos,4,0.5,NORM((fabs)(data.sig[DATA::Don].weight[Signal::Forward]),20.0), (data.sig[DATA::Don].weight[Signal::Forward] >0) ? 5 : 6 );
            if (data.sig[DATA::Acc].weight[Signal::Reverse] != 0 )
              PlotBarF(pos,-4,0.5,NORM((fabs)(data.sig[DATA::Acc].weight[Signal::Reverse]),20.0), (data.sig[DATA::Acc].weight[Signal::Reverse] >0) ? 4 : 7 );
            if (data.sig[DATA::Don].weight[Signal::Reverse] != 0 )
              PlotBarF(pos,-4,0.5,NORM((fabs)(data.sig[DATA::Don].weight[Signal::Reverse]),20.0),(data.sig[DATA::Don].weight[Signal::Reverse] >0) ? 5 : 6 );
        */
    }

}

// ------------------
//  Post analyse
// ------------------
void SensorSignalWAM :: PostAnalyse(Prediction *pred, FILE *MINFO)
{
}
