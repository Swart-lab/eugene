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
#include "Output.h"
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

  int    j, k;
  int    Data_Len = TheSeq->SeqLen;
  DATA	 Data;
  
  DAG Dag(0, Data_Len, PAR,TheSeq);
  
  Dag.LoadDistLength();
  Dag.WeightThePrior();
  
  // --------------------------------------------------------------------------
  // Demarrage de la programmation dynamique
  // --------------------------------------------------------------------------
  for (int nuc = 0; nuc <= Data_Len; nuc++) {
    
    // recuperation des infos
    MSensor->GetInfoAt(TheSeq, nuc, &Data);
    Dag.ShortestPathAlgoForward(nuc,Data);
  }
  
  Dag.WeightThePrior();
  Dag.BuildPrediction();

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
  DNASeq *TheSeq;
  int    Data_Len;
  Prediction *pred;
  char   grnameFile[FILENAME_MAX+1];
  char   grname[FILENAME_MAX+1];
  int    graph;

  fprintf(stderr,"-------------------------------------");
  fprintf(stderr,"--------------------------------\n");

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
	// Preparation sortie graphique + Scores
	// --------------------------------------------------------------------
	if (graph) {
	  int gto       = PAR.getI("Output.gto");
	  int gfrom     = PAR.getI("Output.gfrom");
	  int glen      = PAR.getI("Output.glen");
      
	  // Récupération du nom du fichier d'origine
	  strcpy(grnameFile, BaseName(PAR.getC("fstname")));
	  if ( rindex(grnameFile, '.') != NULL )
	    *rindex(grnameFile, '.') = 0;     // on enleve l'extension
		  
	  // Construction du nom de sortie (*.png)
	  strcpy(grname, PAR.getC("Output.Prefix"));
	  strcat(grname, grnameFile);
	  if(PAR.getC("grnameArg")[0]!='0') { // -ph sans -g ou -g NO pas d'arg
	    strcat(grname,".");
	    strcat(grname,PAR.getC("grnameArg"));
	  }
	  
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

	pred = Predict(TheSeq, MS);
    
	fprintf(stderr,"Optimal path length = %.4f\n",- pred->OptimalPath);
	
	if (graph)
	  pred->plotPred();
    
	Output(TheSeq, MS, pred, sequence, argc, argv, stdout);
	MS->PostAnalyse(pred);
    
	if (graph) {
	  fprintf(stderr,"\nDumping images (\"%s.---.png\")...", grname);
	  fflush(stderr);
	  ClosePNG();
	  fprintf(stderr, "done\n");
	}
	
	// Free used memory
	delete TheSeq; 
	delete MS;
	delete pred;//->resetPred();
	
      } // fin de traitement de chaque séquence....
  
      fprintf(stderr,"-------------------------------------");
      fprintf(stderr,"--------------------------------\n");
  
      return  0;
    }
}
