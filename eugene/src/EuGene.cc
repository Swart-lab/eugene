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
// $Id: EuGene.cc,v 1.109 2015-07-30 14:19:36 sallet Exp $
// ------------------------------------------------------------------
// File:     EuGene.cc
// Contents: This program finds exons/introns and intergenic regions
//           (including UTRs)
// ------------------------------------------------------------------


#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <cstdio>
#include <cstdlib>

#ifdef __APPLE__
// MacOS-X kludge. cmath undefines these macros. Turn them into inlines
#include <math.h>
inline int (isinf)(double r)
{
    return isinf(r);
}
inline int (isnan)(double r)
{
    return isnan(r);
}
#endif

#include <cmath>
#include <cctype>
#include <climits>
#include <cfloat>
#include <ctime>
#include <cassert>
#include <cerrno>
#include <map>
#include <unistd.h>
#include <string.h>
#ifdef HAVE_STRINGS_H
#include <strings.h>
#endif

#include "Const.h"
#include "System.h"
#include "DAG.h"
#include "AltEst.h"
#include "SensorIF.h"
#include "MSensor.h"
#include "Param.h"
#include "DNASeq.h"
#include "BackP.h"
#include "PenaltyDist.h"
#include "Prediction.h"
#include "Parametrization/ParaOptimization.h"

#ifndef FILENAME_MAX
#define FILENAME_MAX        1024
#endif

// ------------------ Globals --------------------
MasterSensor*    MS;
Parameters       PAR;
ParaOptimization OPTIM;
std::vector<short int> activeTracks;              // Vector containing the numbers of the active tracks
std::vector<bool> isActiveTrack(NbTracks, false); // Vector of bool. True if a track is active
std::map<int, vector<int> > bicodingStates;       // Bicoding state content. Ex: bicodingStates[SnglF1R2] = [SnglF1, SnglF2]


// -------------------------------------------------------------------------
// Initialize the active tracks vector.
// At the moment, just copy the values including in the array proActiveTracks
// or eukActiveTracks in Prediction_cte.h
// TODO Load a file and remove these arrays
// -------------------------------------------------------------------------
void InitActiveTracks(int strand)
{
    char* mode = PAR.getC ("EuGene.mode");

    assert(!strcmp(mode, "Prokaryote") || !strcmp(mode, "Eukaryote") || !strcmp(mode, "Prokaryote2") || !strcmp(mode, "Eukaryote2"));

	// Init active track vectors
	activeTracks.clear();
	for (int i=0; i < NbTracks; i++)
		isActiveTrack[i] = false;

    if (!strcmp(mode, "Prokaryote"))
    {
        for (int  i = 0; i < proActiveTracksNb; i++)
        {
            activeTracks.push_back(proActiveTracks[i]);
            isActiveTrack[proActiveTracks[i]] = true;
        }
    }
    else if (!strcmp(mode, "Eukaryote"))
    {
        for (int  i = 0; i < eukActiveTracksNb; i++)
        {
            activeTracks.push_back(eukActiveTracks[i]);
            isActiveTrack[eukActiveTracks[i]] = true;
        }
    }
    else if (!strcmp(mode, "Prokaryote2")) // "Prokaryote2"
    {	
		assert ( (strand == 1) || (strand == -1) );
		if (strand == 1)
		{
			for (int  i = 0; i < proForwardActiveTracksNb; i++)
        	{
            	activeTracks.push_back(proForwardActiveTracks[i]);
            	isActiveTrack[proForwardActiveTracks[i]] = true;
        	}
		}
		else
		{
			for (int  i = 0; i < proReverseActiveTracksNb; i++)
        	{
            	activeTracks.push_back(proReverseActiveTracks[i]);
            	isActiveTrack[proReverseActiveTracks[i]] = true;
        	}
		}
	}
	else if (!strcmp(mode, "Eukaryote2"))
    {
        assert ( (strand == 1) || (strand == -1) );
        if (strand == 1)
        {
            for (int  i = 0; i < eukForwardActiveTracksNb; i++)
            {
                activeTracks.push_back(eukForwardActiveTracks[i]);
                isActiveTrack[eukForwardActiveTracks[i]] = true;
            }
        }
        else
        {
            for (int  i = 0; i < eukReverseActiveTracksNb; i++)
            {
                activeTracks.push_back(eukReverseActiveTracks[i]);
                isActiveTrack[eukReverseActiveTracks[i]] = true;
            }
        }
    }
}

// -------------------------------------------------------------------------
// Compute Optimal Prediction
// -------------------------------------------------------------------------
Prediction* Predict (DNASeq* TheSeq, int From, int To, MasterSensor* MSensor, int strand=0)
{

    int   j, k;
    int   Data_Len = TheSeq->SeqLen;
    DATA	Data;
    int   Forward =  1;//PAR.getI("Sense");
    int   Dir = (Forward ? 1 : -1);

    int   GCVerbose = PAR.getI("EuGene.VerboseGC");
    int   GCLatency = PAR.getI("EuGene.GCLatency");

    // Active the appropriated tracks according to the eugene mode Prokaryote or Eukaryote
    InitActiveTracks(strand);

    DAG   *Dag;
    Prediction *prediction;

    // DynaProg end at the lastpos + 1 to account for final signals.
    int	FirstNuc = (Forward ? From : To+1);
    int   LastNuc  = (Forward ? To+1 : From);

    // Initialize the central DAG
    // one more position here for prediction of partial elements at extremities.
    Dag = new DAG(FirstNuc-Dir, LastNuc+Dir, PAR,TheSeq);

    Dag->LoadDistLength();
    Dag->WeightThePrior();

    // --------------------------------------------------------------------------
    // Demarrage de la programmation dynamique
    // --------------------------------------------------------------------------
    for (int nuc = FirstNuc; nuc != LastNuc+Dir; nuc += Dir)
    {

        // recuperation des infos
        MSensor->GetInfoAt(TheSeq, nuc, &Data);
        if (Forward)
            Dag->ShortestPathAlgoForward(nuc,Data);
        else
            Dag->ShortestPathAlgoBackward(nuc,Data);

        if (nuc && (nuc % GCLatency == 0)) Dag->MarkAndSweep(nuc,GCVerbose,GCLatency);
    }

    Dag->WeightThePrior();
    Dag->BuildPrediction(From,To,Forward);

    Dag->pred->TrimAndUpdate(TheSeq);
    prediction = Dag->pred;

    delete Dag;

    return prediction;
}

// Same with whole sequence by default
Prediction* Predict (DNASeq* TheSeq, MasterSensor* MSensor)
{
    return Predict(TheSeq,0,TheSeq->SeqLen-1,MSensor);
}
// -------------------------------------------------------------------------
// Compute alternative Predictions based on EST
// -------------------------------------------------------------------------
Prediction* AltPredict (DNASeq* TheSeq, int From, int To, MasterSensor* MSensor,
                        AltEst *AltEstDB, Prediction *optpred, int idx)
{
    DATA	Data;
    int   Forward =  1;//PAR.getI("Sense");
    int   Dir = (Forward ? 1 : -1);
    int   GCVerbose = PAR.getI("EuGene.VerboseGC");
    int   GCLatency = PAR.getI("EuGene.GCLatency");
    DAG* Dag;
    Prediction *prediction;

        // Active the appropriated tracks according to the eugene mode Prokaryote or Eukaryote
    InitActiveTracks(0);
    
    // DynaProg end at the lastpos + 1 to account for final signals.
    int	FirstNuc = (Forward ? From : To+1);
    int   LastNuc  = (Forward ? To+1 : From);

    if (!AltEstDB->voae_AltEst[idx].CompatibleWith(optpred))
    {

        Dag = new DAG(FirstNuc-Dir, LastNuc+Dir, PAR,TheSeq);

        Dag->LoadDistLength();
        Dag->WeightThePrior();

        for (int nuc = FirstNuc; nuc != LastNuc+Dir; nuc += Dir)
        {

            // recuperation des infos
            MSensor->GetInfoAt(TheSeq, nuc, &Data);
            AltEstDB->Penalize(idx,nuc,&Data);

            if (Forward)
                Dag->ShortestPathAlgoForward(nuc,Data);
            else
                Dag->ShortestPathAlgoBackward(nuc,Data);

            if (nuc && (nuc % GCLatency == 0)) Dag->MarkAndSweep(nuc,GCVerbose,GCLatency);
        }
        Dag->WeightThePrior();
        Dag->BuildPrediction(From, To, Forward);

        Dag->pred->TrimAndUpdate(TheSeq);
        prediction = Dag->pred;

        delete Dag;
        return prediction;
    }
    return NULL;

}

// -------------------------------------------------------------------------
// Read a fasta file
// -------------------------------------------------------------------------
DNASeq* ReadSequence (char* sequence_name)
{
    DNASeq *TheSeq;

    fprintf(stderr,"-------------------------------------");
    fprintf(stderr,"--------------------------------\nLoading sequence...");
    fflush(stderr);

    TheSeq = new DNASeq(sequence_name);

    fprintf(stderr,"%s, %d bases read, ",TheSeq->Name, TheSeq->SeqLen);

    fprintf (stderr,"GC Proportion = %.1f%%\n",
             (TheSeq->Markov0[BitG]+TheSeq->Markov0[BitC])*100.0);

    return TheSeq;
}



// -------------------------------------------------------------------------
// Fill bicodingStates map.
// TODO: Improve this method: dont have to be global to EuGene.
// -------------------------------------------------------------------------
void InitBicoding()
{
    vector<int> v;
    v.push_back(SnglF1);
    v.push_back(SnglF2);
    bicodingStates[SnglF1F2] = v;
    v.clear();
    v.push_back(SnglF1);
    v.push_back(SnglF3);
    bicodingStates[SnglF1F3] = v;
    v.clear();
    v.push_back(SnglF2);
    v.push_back(SnglF3);
    bicodingStates[SnglF2F3] = v;
    v.clear();
    v.push_back(SnglR1);
    v.push_back(SnglR2);
    bicodingStates[SnglR1R2] = v;
    v.clear();
    v.push_back(SnglR1);
    v.push_back(SnglR3);
    bicodingStates[SnglR1R3] = v;
    v.clear();
    v.push_back(SnglR2);
    v.push_back(SnglR3);
    bicodingStates[SnglR2R3] = v;
    v.clear();
    v.push_back(SnglF1);
    v.push_back(SnglR1);
    bicodingStates[SnglF1R1] = v;
    v.clear();
    v.push_back(SnglF1);
    v.push_back(SnglR2);
    bicodingStates[SnglF1R2] = v;
    v.clear();
    v.push_back(SnglF1);
    v.push_back(SnglR3);
    bicodingStates[SnglF1R3] = v;
    v.clear();
    v.push_back(SnglF2);
    v.push_back(SnglR1);
    bicodingStates[SnglF2R1] = v;
    v.clear();
    v.push_back(SnglF2);
    v.push_back(SnglR2);
    bicodingStates[SnglF2R2] = v;
    v.clear();
    v.push_back(SnglF2);
    v.push_back(SnglR3);
    bicodingStates[SnglF2R3] = v;
    v.clear();
    v.push_back(SnglF3);
    v.push_back(SnglR1);
    bicodingStates[SnglF3R1] = v;
    v.clear();
    v.push_back(SnglF3);
    v.push_back(SnglR2);
    bicodingStates[SnglF3R2] = v;
    v.clear();
    v.push_back(SnglF3);
    v.push_back(SnglR3);
    bicodingStates[SnglF3R3] = v;
    v.clear();
}


// -------------------------------------------------------------------------
//            MAIN
// -------------------------------------------------------------------------
int main  (int argc, char * argv [])
{
    DNASeq     *TheSeq;
    int        Data_Len;
    Prediction *pred;
    FILE       *MISC_INFO;
    char       prefixName[FILENAME_MAX+1];
    char       grname[FILENAME_MAX+1];
    char       miname[FILENAME_MAX+1];
    int        graph;
	bool       debugAltest = true;
	

    fprintf(stderr,"-------------------------------------"
            "--------------------------------\n");

    // Lecture de la ligne d'arg et du fichier .par
    PAR.initParam(argc, argv);

    // init bicoding state map
    InitBicoding();

    if (PAR.getI("ParaOptimization.Use"))
        OPTIM.ParaOptimize(argc, argv);
    else
    {
        // Objectif : limiter les appels à la MAP
        graph = PAR.getI("Output.graph");

        int sequence;
        for (sequence = optind; sequence < argc ; sequence++)
        {
			time_t t1, t2, t3, t4, t5, t6, t7;
			time(&t1); // start analayse the sequence
	

            PAR.set("fstname", argv[sequence]);

            // --------------------------------------------------------------------
            // Lecture de la sequence
            // --------------------------------------------------------------------
            TheSeq = ReadSequence( PAR.getC("fstname") );
            Data_Len = TheSeq->SeqLen;

            // --------------------------------------------------------------------
            // Calcul des positions de prédiction
            // --------------------------------------------------------------------
            int fromPos = PAR.getI("EuGene.from",0,true);
            int toPos   = PAR.getI("EuGene.to",0,true);

            if (toPos)
                toPos = Min(Data_Len-1,toPos);
            else
                toPos = Data_Len-1;

            // --------------------------------------------------------------------
            // Prefix output file name
            // --------------------------------------------------------------------
            strcpy(prefixName, PAR.getC("Output.Prefix"));
            strcat(prefixName, BaseName(PAR.getC("fstname")));
            if ( rindex(prefixName, '.') != NULL )
            {
                if (!strcmp(rindex(prefixName, '.'), ".fasta") ||
                        !strcmp(rindex(prefixName, '.'), ".fsa")   ||
                        !strcmp(rindex(prefixName, '.'), ".tfa")   ||
                        !strcmp(rindex(prefixName, '.'), ".txt"))
                    *rindex(prefixName, '.') = 0;     // on enleve l'extension
            }
            PAR.set("prefixName", prefixName);

            // --------------------------------------------------------------------
            // Preparation sortie graphique + Scores
            // --------------------------------------------------------------------
            if (graph)
            {
                int gto       = PAR.getI("Output.gto");
                int gfrom     = PAR.getI("Output.gfrom");
                int glen      = PAR.getI("Output.glen");

                // Construction du nom de sortie (*.png)
                strcpy(grname, prefixName);

                if ((gfrom <= 0)|| (gfrom >= Data_Len))
                    gfrom = 1;
                if ((gto <= 0)  || (gto <= gfrom) || (gto > Data_Len))
                    gto = Data_Len;

                if ((PAR.getI("Output.gfrom")!=-1) || PAR.getI("Output.gto")!=-1)
                {
                    sprintf(grname+strlen(grname), ".%d", gfrom);
                    sprintf(grname+strlen(grname), "-%d", gto);
                }

                gfrom--;
                gto--;

                if (glen < 0)
                    glen = ((gto-gfrom+1 <= 6000) ? gto-gfrom+1 : 6000);

                InitPNG(PAR.getI("Output.resx"),   PAR.getI("Output.resy"),
                        PAR.getI("Output.offset"), gfrom, gto,
                        PAR.getI("Output.golap"), glen, grname);
            }

            // --------------------------------------------------------------------
            // Init MasterSensor
            // --------------------------------------------------------------------
            MS = new MasterSensor();
            MS->InitMaster(TheSeq);
			time(&t2); // end of Init the sensor

            // --------------------------------------------------------------------
            // Predict: 1st main prediction
            // --------------------------------------------------------------------
            if (PAR.count("AltEst.reference") > 0)
            {
                fprintf(stderr, "Alt mode\n");
                char reffile[FILENAME_MAX+1];
                strcpy(reffile, PAR.getC("AltEst.reference"));
                pred = new Prediction(reffile, TheSeq);
				
				time(&t3); // end of create the Prediction object from the reference GFF3
				
                //pred->Print();
            }
            else
            {
                
                char* mode = PAR.getC ("EuGene.mode");

                if (!strcmp(mode, "Prokaryote2") || !strcmp(mode, "Eukaryote2") )
                {
                    // predict on the forward strand
                    Prediction* pred2 = Predict(TheSeq, fromPos, toPos, MS,  1);
                    // predict on the reverse strand
                    pred = Predict(TheSeq, fromPos, toPos, MS, -1);
                    pred->AppendPred(pred2);
                }
                else
                {
                    pred = Predict(TheSeq, fromPos, toPos, MS);
                    //pred->Print();
                    fprintf(stderr,"Optimal path length = %.4f\n",- pred->optimalPath);
                }
                //pred->Print();
            }

            // --------------------------------------------------------------------
            // Textual and graphical output
            // --------------------------------------------------------------------
            if (graph)
                pred->PlotPred();

            if ( ! PAR.getI("AltEst.use") )
                pred->Print(TheSeq, MS);

            strcpy(miname, prefixName);
            MISC_INFO = FileOpen(NULL, strcat(miname, ".misc_info"), "wb");
            pred->PrintGeneInfo(MISC_INFO);
            MS->PostAnalyse(pred, MISC_INFO);
            fclose(MISC_INFO);


            if (graph)
            {
                fprintf(stderr,"Dumping images (\"%s.---.png\")...", grname);
                fflush(stderr);
                ClosePNG();
                fprintf(stderr, "done\n");
            }

			time(&t4); // end of MS PostAnalyse
            // --------------------------------------------------------------------
            // Load Alternative EST data (if any)
            // --------------------------------------------------------------------
            if (PAR.getI("AltEst.use"))
            {
                int ExonBorderMatchThreshold = PAR.getI("AltEst.ExonBorderMatchThreshold");
                int RepredictMargin          = PAR.getI("AltEst.RepredictMargin");
	
                int newGene = 0; // if a splice variant has no base gene, it is a "new" gene. counter needed for gene number
                

                AltEst *AltEstDB = new AltEst(TheSeq);
				time(&t5); // end of AltEst build
				
                std::vector <Prediction*> vPred;
                Prediction*               AltPred;
                Gene*                     baseGene;
				
				time_t depart, arrivee;
				time(&depart);
                for (int altidx = 0; altidx < AltEstDB->totalAltEstNumber; altidx++)
                {
                    int localFrom,localTo;
					if (debugAltest && altidx%10000 == 0)
					{
						if (altidx > 0) 
						{
							time(&arrivee);
							double diff = difftime(arrivee, depart);
							fprintf(stderr, "[ALT EST] Time to perform 10000 AltPredict (index=%d) %.f seconds (%d minutes)\n", altidx, diff, int(diff/60));
						}
						time(&depart);
					}
					
                    localFrom = Max(fromPos, AltEstDB->voae_AltEst[altidx].GetStart()-RepredictMargin);
                    localTo   = Min(toPos,   AltEstDB->voae_AltEst[altidx].GetEnd()+RepredictMargin);

                    AltPred = AltPredict(TheSeq,localFrom,localTo,MS,AltEstDB,pred,altidx);

                    if (AltPred)
                    {
                        if ( (AltPred->vGene[0]->cdsStart == -1) || (AltPred->vGene[0]->cdsEnd == -1))
                        {
                            delete AltPred;
                            continue;
                        }
                        // Delete the gene of the alt prediction which doesn't overlap the EST
                        AltPred->DeleteOutOfRange(AltEstDB->voae_AltEst[altidx].GetStart(),AltEstDB->voae_AltEst[altidx].GetEnd(), AltEstDB->voae_AltEst[altidx].GetStrand());
                        // If genes overlapping the EST was found and if the prediction is original
                        if ( (AltPred->nbGene > 0) && (AltPred->IsOriginal(pred,vPred,ExonBorderMatchThreshold)) )
                        {
                            fprintf(stderr,"Optimal path length = %.4f\n",- AltPred->optimalPath);
                            baseGene = pred->FindGene(AltPred->vGene[0]->trStart,AltPred->vGene[0]->trEnd, AltPred->vGene[0]->GetStrand());
                            if (baseGene)
                            {
                                baseGene->hasvariant++;
                                AltPred->vGene[0]->isvariant = true;
                                AltPred->vGene[0]->hasvariant = baseGene->hasvariant;
                                AltPred->vGene[0]->geneNumber = baseGene->geneNumber;
                                baseGene->tuStart = ( baseGene->tuStart ) ? Min(baseGene->tuStart,AltPred->vGene[0]->trStart)
                                                    : Min(baseGene->trStart,AltPred->vGene[0]->trStart);
                                baseGene->tuEnd   = ( baseGene->tuEnd )   ? Max(baseGene->tuEnd,AltPred->vGene[0]->trEnd)
                                                    : Max(baseGene->trEnd,AltPred->vGene[0]->trEnd);
                            }
                            else
                            {
                                fprintf(stderr,"New gene predicted by alternative spliced gene prediction.\n");
                                AltPred->vGene[0]->geneNumber = pred->nbGene + newGene++;
                            }
                            vPred.push_back(AltPred);
                        }
                        else delete AltPred;
                    }
                } 
                
                time(&t6); // end of analyse of all altEst
                pred->Print(TheSeq, MS);
                for (int idx = 0; idx < vPred.size(); idx++)
                {
                    vPred[idx]->Print(TheSeq, MS,NULL,1);
                }
                time(&t7); // end of print all the predictions
                //delete AltPred;
				if (debugAltest)
				{
					double part1= difftime(t2, t1);
					fprintf(stderr, "[ALT EST] Duration of start until the init sensor: %.f seconds (%d minutes)\n", part1, int (part1/60));
					double part2 = difftime(t3, t2);
					fprintf(stderr, "[ALT EST] Duration of Prediction object creation from the reference GFF3: %.f seconds (%d minutes)\n", part2, int (part2/60));
					double part3 = difftime(t4, t3);
					fprintf(stderr, "[ALT EST] Duration of the MasterSensor Post analyse: %.f seconds (%d minutes)\n", part3, int (part3/60));
					double part4= difftime(t5, t4);
					fprintf(stderr, "[ALT EST] Duration of the AltEstDB building: %.f seconds (%d minutes)\n", part4, int (part4/60));
					fprintf(stderr, "[ALT EST] Keep %d AltEst\n", AltEstDB->totalAltEstNumber);
					double part5= difftime(t6, t5);
					fprintf(stderr, "[ALT EST] Duration of AltPredict:  %.f seconds (%d minutes)\n", part5, int(part5/60));
					double part6= difftime(t7, t6);
					fprintf(stderr, "[ALT EST] Duration of printing of all predictions: %.f seconds (%d minutes)\n", part6, int(part6/60));
					double allparts = difftime(t7, t1);
					fprintf(stderr, "[ALT EST] Total duration : %.f seconds (%d minutes)\n", allparts, int(allparts/60));
				}
            }

            // Free used memory
            delete TheSeq;
            delete MS;
            delete pred;

            fflush(stderr);
            fflush(stdout);
        } // fin de traitement de chaque séquence....

        fprintf(stderr,"-------------------------------------"	      "--------------------------------\n");

        return  0;
    }
	
}
