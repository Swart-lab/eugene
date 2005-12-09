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
// File:     EuGene.cc
// Contents: This program finds exons/introns and intergenic regions 
//           (including UTR)
// ------------------------------------------------------------------


#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <cstdio>
#include <cstdlib>

#ifdef __APPLE__
// MacOS-X kludge. cmath undefines these macros. Turn them into inlines 
#include <math.h>
inline int (isinf)(double r) { return isinf(r); }
inline int (isnan)(double r) { return isnan(r); }
#endif

#include <cmath>
#include <cctype>
#include <climits>
#include <cfloat>
#include <ctime>
#include <cassert>
#include <cerrno>
#include <unistd.h>
#include <string.h>
#ifdef HAVE_STRINGS_H
#include <strings.h>
#endif

#include "Const.h"
#include "System.h"
#include "DAG.h"
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

// -------------------------------------------------------------------------
// Predict genes
// -------------------------------------------------------------------------
Prediction* Predict (DNASeq* TheSeq, MasterSensor* MSensor)
{

  int   j, k;
  int   Data_Len = TheSeq->SeqLen;
  DATA	Data;
  int   Forward =  1;//PAR.getI("Sense"); 
  int   Dir = (Forward ? 1 : -1);
  int	FirstNuc = (Forward ? 0 : Data_Len);
  int   LastNuc  = (Forward ? Data_Len : 0)+Dir;
  int   GCVerbose = PAR.getI("EuGene.VerboseGC");
  int   GCLatency = PAR.getI("EuGene.GCLatency");
  DAG   Dag(FirstNuc-Dir, LastNuc, PAR,TheSeq);
  
  Dag.LoadDistLength();
  Dag.WeightThePrior();

  // --------------------------------------------------------------------------
  // Demarrage de la programmation dynamique
  // --------------------------------------------------------------------------
  for (int nuc = FirstNuc; nuc != LastNuc; nuc += Dir) {
    
    // recuperation des infos
    MSensor->GetInfoAt(TheSeq, nuc, &Data);
    if (Forward) 
      Dag.ShortestPathAlgoForward(nuc,Data);
    else
      Dag.ShortestPathAlgoBackward(nuc,Data);

    if (nuc && (nuc % GCLatency == 0)) Dag.MarkAndSweep(nuc,GCVerbose);
  }

  Dag.WeightThePrior();
  Dag.BuildPrediction(Forward);

  return Dag.pred;
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

  fprintf(stderr,"-------------------------------------"
	  "--------------------------------\n");

  // Lecture de la ligne d'arg et du fichier .par
  PAR.initParam(argc, argv);

  if (PAR.getI("ParaOptimization.Use")) 
    OPTIM.ParaOptimize(argc, argv);
  else 
    {
      // Objectif : limiter les appels à la MAP
      graph = PAR.getI("Output.graph");

      int sequence;
      for (sequence = optind; sequence < argc ; sequence++) {
	
	PAR.set("fstname", argv[sequence]);

	// --------------------------------------------------------------------
	// Lecture de la sequence    
	// --------------------------------------------------------------------
	TheSeq = ReadSequence( PAR.getC("fstname") );
	Data_Len = TheSeq->SeqLen;
    
	// --------------------------------------------------------------------
	// Prefix output file name
	// --------------------------------------------------------------------
	strcpy(prefixName, PAR.getC("Output.Prefix"));
	strcat(prefixName, BaseName(PAR.getC("fstname")));
	if ( rindex(prefixName, '.') != NULL ) {
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
	if (graph) {
	  int gto       = PAR.getI("Output.gto");
	  int gfrom     = PAR.getI("Output.gfrom");
	  int glen      = PAR.getI("Output.glen");
      		  
	  // Construction du nom de sortie (*.png)
	  strcpy(grname, prefixName);
	  
	  if ((gfrom <= 0)|| (gfrom >= Data_Len))
	    gfrom = 1;
	  if ((gto <= 0)  || (gto <= gfrom) || (gto > Data_Len))
	    gto = Data_Len;
	  
	  if ((PAR.getI("Output.gfrom")!=-1) || PAR.getI("Output.gto")!=-1) {
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

	// --------------------------------------------------------------------
	// Predict
	// --------------------------------------------------------------------
	pred = Predict(TheSeq, MS);
       	fprintf(stderr,"Optimal path length = %.4f\n",- pred->optimalPath);
	
	// --------------------------------------------------------------------
	// Textual and graphical output
	// --------------------------------------------------------------------
	if (graph)
	  pred->PlotPred();

	pred->Print(TheSeq, MS);

	strcpy(miname, prefixName);
	MISC_INFO = FileOpen(NULL, strcat(miname, ".misc_info"), "wb");
	pred->PrintGeneInfo(MISC_INFO);
	MS->PostAnalyse(pred, MISC_INFO);
	fclose(MISC_INFO);

    	if (graph) {
	  fprintf(stderr,"Dumping images (\"%s.---.png\")...", grname);
	  fflush(stderr);
	  ClosePNG();
	  fprintf(stderr, "done\n");
	}

	// Free used memory
	delete TheSeq; 
	delete MS;
	delete pred;
	
	fflush(stderr);
	fflush(stdout);
      } // fin de traitement de chaque séquence....
	    
      fprintf(stderr,"-------------------------------------"
	      "--------------------------------\n");
    
      return  0;
    }
}
