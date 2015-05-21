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
# include <dirent.h>
# include <sys/types.h>
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
    char* signalType = PAR.getC("SignalWAM.type",          GetNumber()); // start, acceptor or donor
    char* species    = PAR.getC("SignalWAM.species",       GetNumber()); // species dir models/WAM/species has to exist
    scoringMode      = PAR.getI("SignalWAM.scoringMode",   GetNumber());
    
    /////////////////////////////////////////////////////////////////////////////////////
    // Assert signalType = donor, acceptor or start
    /////////////////////////////////////////////////////////////////////////////////////
    
    int signalLen;   
    if (strcmp(signalType, "donor") == 0 || strcmp(signalType, "acceptor") == 0)
    {
        signalLen = 2;
    }
    else if (strcmp(signalType, "start") == 0)
    {
        // start
        signalLen = 3;
    }
    else
    {
        fprintf(stderr, "Error \"%s\" undefined type in the parameter file"
        " for SignalWAM.type[%d]\n", signalType, GetNumber());
        fprintf(stderr, "Type must be : start, acceptor or donor.\n");
        exit(2);  
    }
    
    /////////////////////////////////////////////////////////////////////////////////////
    // In the wam directory  $EUGENEDIR/models/wam/species/signalType
    // Search the best available wam models (max order)
    /////////////////////////////////////////////////////////////////////////////////////
    
    char *wamdir = new char[FILENAME_MAX+1];
    sprintf(wamdir, "%s/%s/%s/%s/", PAR.getC("eugene_dir"), WAM_DIR, species, signalType); 
    DIR *dir = opendir(wamdir);
    if (!dir)
    {
        fprintf(stderr, "Cannot open wam model directory %s\n", wamdir);
        exit(2);
    }
    
    // Search wam models of highest order:
    // Sort files in descending, then get the first file name (file of highest order) 
    struct dirent *entry = NULL;
    std::vector<string> vFiles;
    while ( (entry = readdir(dir)) != NULL )
    {
        if (strcmp(entry->d_name, ".") != 0 && strcmp(entry->d_name, "..") != 0)
        {           
            vFiles.push_back(entry->d_name);
        }
    }
    if (vFiles.size() <= 0)
    {
        fprintf(stderr, "No wam models found.\n");
        exit(2); 
    }
    sort(vFiles.rbegin(), vFiles.rend());
    string wamFileName = vFiles.at(0); // WAM.o02.uc020.dc020.AG.TP.39
    
    
    ///////////////////////////////////////////////////////////////////////////////////////
    // Extract information from the file named WAM.o02.uc020.dc020.AG.TP.39 :
    // wam order, upstream and downstream lengths  and site nucleotides
    /////////////////////////////////////////////////////////////////////////////////////
      
    int markovianOrder = atoi(wamFileName.substr(5,2).c_str());  //  02
    upstreamLen        = atoi(wamFileName.substr(10,3).c_str()); // 020
    downstreamLen      = atoi(wamFileName.substr(16,3).c_str()); // 020
    pattern            = new char[4];
    strcpy(pattern, wamFileName.substr(20,signalLen).c_str());   // AG
    for(int i=0; i<strlen(pattern); i++)
        pattern[i] = tolower(pattern[i]);
    
    char *modelfilename = new char[FILENAME_MAX+1];
    strcat(modelfilename, wamdir);
    strcat(modelfilename, wamFileName.substr(0,20+signalLen).c_str());
    //printf ("Fichier: %s =>Order: %i, uc %i, dc %i, patt %s, model file %s\n", wamFileName.c_str(), 
    //           markovianOrder,upstreamLen, downstreamLen, pattern, modelfilename ) ;

    ///////////////////////////////////////////////////////////////////////////////////////
    // Initialization according to the signaltype and the scoringMode
    // codingside = 1: score only the coding side around the site
    /////////////////////////////////////////////////////////////////////////////////////
    
    motifLength       = upstreamLen + signalLen + downstreamLen;    
    scoringMotifStart = 1;
    scoringMotifEnd   = motifLength;

    if (!strcmp(signalType, "start"))
    {
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
        sigType     = DATA::Acc;
        type        = Type_Acc;
        newStatePos = 3;
        if (scoringMode == 1) //codingside
        {
            scoringMotifStart = upstreamLen+3;
            scoringMotifEnd   = motifLength;
        }
    }
    
    ///////////////////////////////////////////////////////////////////////////////////////
    // Build WAM object and launch signal search
    /////////////////////////////////////////////////////////////////////////////////////
    
    WAModel = new WAM(markovianOrder, motifLength, "ACGT", modelfilename);
    
    SearchSignal(X);
    
    
    delete[] wamdir;
    delete[] modelfilename;
    delete dir;
    
}

// ----------------------
//  Default destructor.
// ----------------------
SensorSignalWAM :: ~SensorSignalWAM()
{
    delete [] WAModel;
    delete [] pattern;
    vPosF.clear();
    vPosR.clear();
    vScoreF.clear();
    vScoreR.clear();
}

// ---------------------
//  Init
// ----------------------
void SensorSignalWAM :: Init(DNASeq *X)
{
    scaleCoef    = PAR.getD("SignalWAM.scaleCoef*", GetNumber());
    scalePenalty = PAR.getD("SignalWAM.scalePenalty*", GetNumber());

    positionGiveInfo = -1;
    indexR = indexF = 0;
    
    if (PAR.getI("Output.graph"))
        Plot(X);
}



void SensorSignalWAM :: SearchSignal(DNASeq *X)
{
    int i, j;
    char* Site = new char[motifLength + 1];
    Site[motifLength] = '\0';
    
    int sigStartPos;
    int sigStopPos;
    bool isSigF, isSigR;
    int sigLength = strlen(pattern);
    
    for (int pos = 0; pos <= X->SeqLen; pos++)
    {
        //  ------------------------------------------------------------
        // Forward strand
        // -----------------------------------------------------------------
        isSigF      = true;
        sigStartPos = pos - newStatePos + 1;
        sigStopPos  = sigStartPos + sigLength - 1;
        
        // check there is a signal at the current position on forward strand
        for (i = 0; i < sigLength; i++)
        {
            if ( (*X)[sigStartPos + i] != pattern[i] )
            {
                isSigF = false;
                break;
            }
        }
        
        // If signal detects, extract the context sequence and compute the score
        if ((isSigF) && (sigStartPos - upstreamLen > 0) && (sigStopPos + downstreamLen < X->SeqLen))
        {
            j = 0;
            // extract motif around the signal
            for (i = sigStartPos - upstreamLen; i <= sigStopPos + downstreamLen; i++)
            {
                Site[j] = (*X)[i];
                j++;
            }

            vPosF.push_back(pos);
            vScoreF.push_back(WAModel->ScoreTheMotif(Site, scoringMotifStart, scoringMotifEnd));
        }
        
        // -----------------------------------------------------------------
        // Reverse strand
        // ------------------------------------------------------------------
        isSigR      = true;
        sigStartPos = pos + newStatePos - 2;
        sigStopPos  = sigStartPos - sigLength + 1;
        
        // check there is a signal at the current position on reverse strand
        for (i = sigLength - 1; i > -1; i--)
        {
            if ( (*X)(sigStartPos - i) != pattern[i] )
            {
                isSigR = false;
                break;
            }
        }

        
        // If signal detects, extract the context sequence and compute the score
        if ((isSigR) && (sigStartPos + upstreamLen < X->SeqLen) && (sigStopPos - downstreamLen > 0))
        {
            j = 0;           
            for (i = sigStartPos + upstreamLen; i >= sigStopPos - downstreamLen; i--)
            {
                Site[j] = (*X)(i);
                j++;
            }       
            //if (scoringMode == 2) score = WAModel->ScoreTheSubMotif(Site, 1, upstreamLen) * WAModel->ScoreTheSubMotif(Site, upstreamLen+sigLength+1, motifLength);
            
            vPosR.push_back(pos);
            vScoreR.push_back(WAModel->ScoreTheMotif(Site, scoringMotifStart, scoringMotifEnd));
        }
    }
    
    delete [] Site;     
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
    bool update = false;
    if ( (positionGiveInfo == -1) || (pos != positionGiveInfo+1) ) update = true; // update indexes on vectors
    positionGiveInfo = pos;
    
    if ( !vPosF.empty() ) 
    {
        if (update) 
            indexF = lower_bound(vPosF.begin(), vPosF.end(), pos) - vPosF.begin();
        
        if ( indexF < (int)vPosF.size() && vPosF[indexF] == pos ) 
        {
            d->sig[sigType].weight[Signal::Forward] += (scaleCoef * vScoreF[indexF]) + scalePenalty;
            indexF++;
        }
    }
    
    if ( !vPosR.empty() ) 
    {
        if (update) 
            indexR = lower_bound(vPosR.begin(), vPosR.end(), pos) - vPosR.begin();
        
        if( indexR < (int)vPosR.size() && vPosR[indexR] == pos ) 
        {
            d->sig[sigType].weight[Signal::Reverse] += (scaleCoef * vScoreR[indexR]) + scalePenalty;
            indexR++;
        }
    }
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
