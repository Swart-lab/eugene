// This program build a Weight Array Model (WAM) for a WAM EuGene sensor.
// The goal is to model a signal in a genomic sequence (e.g. splice sites).
// Build the model from a file containing short signal sequences (multi-fasta format).
// Stock results in L (L=WAM length) binary files : one per position in the signal.
// g++ WAMbuilder.cc -o WAMbuild

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "../../0_SensorTk/markov.h"
#include "../../0_SensorTk/markov.cc"

//------ REQUIRED OPTIONS ----------------//
//
//int WAMorder  = 1;
//
//int WAMlength = 9;  // always < 100 !!
//
//int consbeg = 2; // the consensus begins at nt n�consbeg
//
//char* signalfile   = "ARA.DON.12.fasta";
//
//char* outname     = "WAM.ARA.DON.L9.TP.";
//
//----------------------------------------//

//char* signalfile   = "ARA.START9.FP.fasta";
//char* outname     = "WAM.ARA.START9.FP.";

//int WAMlength = 7;  // always < 100 !!
//int consbeg = 4; // the consensus begins at nt n�consbeg
//char* signalfile   = "ARA.ACC.12.fasta";
//char* outname     = "WAM.ARA.ACC.L7.TP.";
//char* signalfile   = "ARA.ACC.12.FP.fasta";
//char* outname     = "WAM.ARA.ACC.L7.FP.";

int main(int argc, char *argv[])
{
  int i,j,k;
  FILE* file;

  if (argc != 6) {
    fprintf(stderr," 
This program build a Weight Array Model (WAM) for a WAM EuGene sensor.
The goal is to model a signal in a genomic sequence (e.g. splice sites).
Build the model from a file containing short signal sequences (multi-fasta format).
Stock results in L (L=WAM length) binary files : one per position in the signal.
\nusage (5 arguments) : %s WAMorder  WAMlength  ConsensusStart  signalinputfile  outname\nexemple: %s 1 9 2 ARA.DON.12.fasta WAM.ARA.DON.L9.TP.\n",argv[0],argv[0]);
    exit(1);
  }
  int WAMorder= atoi(argv[1]);
  int WAMlength = atoi(argv[2]);
  int consbeg = atoi(argv[3]);
  char* signalfile = argv[4];
  char* outname  = argv[5];

  int basenamelen= strlen(outname);

  ChaineADN* ADN = new ChaineADN;

  TabChaine<ChaineADN,int>* WAMCOUNT;
  TabChaine<ChaineADN,unsigned short int>* WAMMOD;

  WAMCOUNT = new TabChaine<ChaineADN,int>[WAMlength](WAMorder,ADN);
  WAMMOD   = new TabChaine<ChaineADN,unsigned short int>[WAMlength](WAMorder,ADN);

  //----------------------------------------------------------------//
  // COUNT words occurences for each position
  fprintf(stderr," - counting words occurences (order %d)  ... position...",WAMorder);
  for (i=0;i<WAMlength;i++){
    file= fopen (signalfile,"rt");
    if (file==NULL) {
      printf("cannot open consensus sequences file %s\n",signalfile);
      exit(1);
    }
    fprintf(stderr,"%d...",i);
    WAMCOUNT[i].fichier2compte(file,consbeg+i-1,consbeg+i);
    fclose(file);
  }
  fprintf(stderr,"done\n");

//  for (i=0;i<WAMlength;i++){
//    printf("position %d:\n",i);
//    WAMCOUNT[i].affichage(0);
//  }

  //----------------------------------------------------------------//
  // ADD pseudocount to each word occurence (to avoid null probs)
  fprintf(stderr," - adding pseudocounts...");
  for (i=0;i<WAMlength;i++){
    WAMCOUNT[i].pseudocount();
  }
  fprintf(stderr,"done\n");
  //  WAMCOUNT.affichage(0);

  //----------------------------------------------------------------//
  // COMPUTE PROBS (store results in matrix WAMMOD)
  fprintf(stderr," - computing probs from frequences for each position...");
  for (i=0;i<WAMlength;i++){
    WAMMOD[i].compte2probas(WAMCOUNT[i]);
    fprintf(stderr,"%d...",i);
  }
  fprintf(stderr,"done\n");

  //----------------------------------------------------------------//
  // SAVE data in the matrix file
  fprintf(stderr," - saving data in matrix files... %s...",outname);

  for (i=0;i<WAMlength;i++) {
    char* name = new char[(basenamelen+3)];
    name[basenamelen+3]='\0';
    sprintf(name,"%s",outname);
    if (i<10) 
      sprintf(name+basenamelen,"0%d",i);
    else  
      sprintf(name+basenamelen,"%d",i);

    file= fopen (name,"wb");
    if (file==NULL) {
      printf("cannot open matrix file %s\n",name);
      exit(1);
    }
    WAMMOD[i].sauve2fichier(file);
    fclose(file);
    fprintf(stderr,"%d...",i);
    delete name;
  }
  fprintf(stderr,"done\n");

  //----------------------------------------------------------------//
  // PRINT data (stdout)
  printf("----- WAM COUNT -----\n");
  for (i=0;i<WAMlength;i++) {
    printf("WAMCOUNT[%d]:\n",i);
    WAMCOUNT[i].affichagevaleurs();
  }

/*
  //----------------------------------------------------------------//
  // LOAD data from a matrix file
  //  WAMMOD.initialisation();
  //file= fopen ("../Markov/sp.modele.dou.o2.bin","rb");
  //if (file==NULL) printf("cannot open matrix file\n");
  //WAMMOD.chargefichier(file);
  //fclose (file);
  //WAMMOD.affichage(0);
  // g++ WAMbuilder.cc -o WAMbuild
*/  
for (i=0;i<WAMlength;i++){
  WAMCOUNT[i].~TabChaine();
  WAMMOD[i].~TabChaine();
}
delete ADN;
}
