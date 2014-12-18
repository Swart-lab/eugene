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
// File:     WAM.cc
// Contents: class WAM
// ------------------------------------------------------------------

#include <ctype.h>

#include "WAM.h"

#include "markov.cc"

/*******************************************************
 **                  WAM                   **
 *******************************************************/

// ----------------------
//  Default constructor
// ----------------------
WAM :: WAM ()
{
  MarkovianOrder = 0;
  MotifLength    = 0;
  Alphabet       = new Chaine("ACGT");
  
  for (int i=0; i< MotifLength-MarkovianOrder; i++) {
    TPMOD.push_back( new TabChaine<Chaine,unsigned short int>(MarkovianOrder,Alphabet) );
    FPMOD.push_back( new TabChaine<Chaine,unsigned short int>(MarkovianOrder,Alphabet) );
  }
}

// --------------------------------------------//
//           Constructor
// --------------------------------------------//
WAM :: WAM (int order, int length, const char* alphabet, char* prefixfilename)
{
  int i;
  int prefixnamelength;
  char* filename;
  FILE* fp;
  char TPfile[FILENAME_MAX+1];
  char FPfile[FILENAME_MAX+1];
  MarkovianOrder = order;
  MotifLength    = length;
  Alphabet       = new Chaine(alphabet);
  
  for (int i=0; i<MotifLength-MarkovianOrder; i++) {
    TPMOD.push_back( new TabChaine<Chaine,unsigned short int>(MarkovianOrder,Alphabet) );
    FPMOD.push_back( new TabChaine<Chaine,unsigned short int>(MarkovianOrder,Alphabet) );
  }

  prefixnamelength = strlen(prefixfilename);
  strcpy(TPfile,prefixfilename);
  strcpy(FPfile,prefixfilename);
  strcat(TPfile,TPFILESUFFIX);
  strcat(FPfile,FPFILESUFFIX);

// Loading models
// (=binary matrix files containing a markov model, one per position of each motif)
  fprintf (stderr,"Reading WAM models...  ");
  fflush (stderr);
  for (i=0; i< MotifLength-MarkovianOrder ;i++) {
    fprintf (stderr,"%d ",i);
    fflush (stderr);
    filename= new char[FILENAME_MAX+1];
    sprintf(filename,"%s",TPfile);
    if (i<10) sprintf(filename+prefixnamelength+SUFFIXLENGTH,"0%d",i);
    else sprintf(filename+prefixnamelength+SUFFIXLENGTH,"%d",i);
    fp=fopen(filename,"rb");
    if  (!fp) {
      fprintf (stderr, "ERROR:  in WAM.cc : could not open file %s \n", filename);
      exit (1);
    }
    if (TPMOD[i]->chargefichier(fp)) {
      fprintf(stderr,"Error when reading model file %s\n",filename);
      exit(2);
    }
    fclose (fp);
    delete[] filename;
    
    filename= new char[FILENAME_MAX+1];
    sprintf(filename,"%s",FPfile);
    if (i<10) sprintf(filename+prefixnamelength+SUFFIXLENGTH,"0%d",i);
    else sprintf(filename+prefixnamelength+SUFFIXLENGTH,"%d",i);
    fp=fopen(filename,"rb");
    if  (!fp) {
      fprintf (stderr, "ERROR:  in WAM.cc : could not open file %s \n", filename);
      exit (1);
    }
    if (FPMOD[i]->chargefichier(fp)) {
      fprintf(stderr,"Error when reading model file %s\n",filename);
      exit(2);
    }
    fclose (fp);
    delete[] filename;
  }
  fprintf (stderr,"... done\n");
  
  /*
  fprintf(stdout, "TP models:\n");
   for (i=0;i<MotifLength-MarkovianOrder;i++){
     fprintf(stdout, "position %d:\n",i);
     TPMOD[i]->affichage(0);
   }
  fprintf(stdout, "FP models:\n");
   for (i=0;i<MotifLength-MarkovianOrder;i++){
     fprintf(stdout, "position %d:\n",i);
     FPMOD[i]->affichage(0);
   }
   */
     
  fflush (stderr);
}

// ----------------------
//  Default destructor.
// ----------------------
WAM :: ~WAM ()
{
  delete Alphabet;
  for (unsigned int i=0; i<TPMOD.size(); i++) {
    delete TPMOD[i];
    delete FPMOD[i];
  }
}

double WAM :: ScoreTheMotif (char* motif) 
// motif must include an amount context of MarkovianOrder length
{
  bool debug = false;
  int i, j;
  int outofalphabet=0;
  double score=0.0;
  if (debug) fprintf(stdout, "motif %s\n", motif);
  
  char* word = new char[MarkovianOrder+2];
  word[MarkovianOrder+1] = '\0'; // word stores a short word for asking markovian probs

  //  at each position of the motif
  for (i=0 ; i< MotifLength-MarkovianOrder ; i++) {
    for (j=0;j<=MarkovianOrder;j++) { // Read a word of MarkovianOrder length 
      word[j] = toupper(motif[i+j]);
      // test if an unknown letter is present, out of the alphabet (like "n" for nucleotid)
      if ((unsigned)Alphabet->operator[](word[j]) == Alphabet->taille) {
	outofalphabet=1;
      }
    }
    if (debug) fprintf(stdout, "--> i=%d, word=%s\t", i, word);
    if (outofalphabet == 0) {
      // likelihood ratio: log ( proba(nt with true model)/proba(nt with false model) )
        if (debug) fprintf(stdout, "score+= (%f - %f)", TPMOD[i]->usi2real(TPMOD[i]->proba(word,MarkovianOrder)), FPMOD[i]->usi2real(FPMOD[i]->proba(word,MarkovianOrder)));
      score += 
	log (TPMOD[i]->usi2real(TPMOD[i]->proba(word,MarkovianOrder)))  -
	log (FPMOD[i]->usi2real(FPMOD[i]->proba(word,MarkovianOrder)))  ;
    }
    if (debug) fprintf(stdout, "\n");
    outofalphabet=0; 
  }

  delete []  word;
  return score;
}
