// This program can build a proteic markov model file for EuGeneAS
// from a proteic data base file in multifasta format (typicaly SwissProt).

// A proteic markov model is generaly less accurate than a nucleic one,
// but is not species-specific, so it can be used as a general coding model,
// without requiring a species-specific genes training data set.
// A set of species specific proteins can also be used.
// For more explanation, see Foissac's PhD thesis (available in 2004 I hope...).

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "markov.h"
#include "System.h"
#include "System.cc"
#include "DNASeq.h"
#include "DNASeq.cc"

int ORDREMAX = 4;  /* maximum of 4 recommended ( pow(21,5) > 4 000 000 ...) */
char* DBfile     = "/www/db/sp";                 // infile
char* MATRIXfile = "swissprot.ordre4.usi.bin";  // outfile
//char* DBfile     = "/Ange/sylvain/swissprot.arath";
//char* MATRIXfile = "swissprot.arath.ordre4.usi.bin";

int main()
{
  int i,j,k;
  FILE* file;
  TabChaine<ChainePROT21,int> COUNT (ORDREMAX,new ChainePROT21);
  COUNT.initialisation();
  TabChaine<ChainePROT21,unsigned short int> MOD(ORDREMAX,new ChainePROT21);
  MOD.initialisation();
  //  TabChaine<ChainePROT21,double> MOD (ORDREMAX,new ChainePROT21);

  //----------------------------------------------------------------//
  // READ the PROTEIC DBfile file (multifasta, e.g. SwissProt)
  fprintf(stderr," - opening proteic DB file...");
  file= fopen (DBfile,"rt");
  if (file==NULL) {
    printf("cannot open swissprot file %s\n",DBfile);
    exit(1);
  }
  fprintf(stderr,"done\n");

  //----------------------------------------------------------------//
  // COUNT words occurences
  fprintf(stderr," - counting words occurences...");
  COUNT.fichier2compte(file);
  fclose(file);
  fprintf(stderr,"done\n");

  //----------------------------------------------------------------//
  // ADD 1 pseudocount to each word occurence (to avoid null probs)
  fprintf(stderr," - adding pseudocounts...");
  COUNT.pseudocount(1);
  fprintf(stderr,"done\n");
  //  COUNT.affichage(0);

  //----------------------------------------------------------------//
  // COMPUTE PROBS (store results in matrix MOD)
  fprintf(stderr," - computing probs...");
  MOD.compte2probas(COUNT);
  fprintf(stderr,"done\n");
  //  MOD.affichage(0);

  //----------------------------------------------------------------//
  // SAVE data in the matrix file
  fprintf(stderr," - saving data in matrix file...");
  file= fopen (MATRIXfile,"wb");
  if (file==NULL) {
    printf("cannot open sp matrix file %s\n",MATRIXfile);
    exit(1);
  }
  MOD.sauve2fichier(file);
  fclose(file);
  fprintf(stderr,"done\n");

  //----------------------------------------------------------------//
  // LOAD data from a matrix file
  //  MOD.initialisation();
  //file= fopen ("../Markov/sp.modele.dou.o2.bin","rb");
  //if (file==NULL) printf("cannot open matrix file\n");
  //MOD.chargefichier(file);
  //fclose (file);
  //MOD.affichage(0);
  // g++ buildMarkovProt.cc -o buildMarkovProt
}
