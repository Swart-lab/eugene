#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <string.h>
#ifdef HAVE_STRINGS_H
#include <strings.h>
#endif

#include "MSensor.h"
#include "Param.h"
#include "DNASeq.h"
#include "Hits.h"
#include "Const.h"

void Output (DNASeq *, char *Choice, int sequence, int argc, char * argv[]);

int HitsCompareLex(const void *, const void *);

int IsPhaseOn(char, int);


// Convertit les phases 0-6 en 1 2 3 -1 -2 -3 0
int PhaseAdapt(char);

void PrintPhase(char);

// Dump les signaux au format util user
void DumpSignals(int Len, REAL** Donor, REAL **Acceptor,
		 REAL** Stop, REAL **Start,FILE* flot);


// Verif coherence EST: calcul le nombre de nuc. coherents et
// incoherents avec les match est
// debut/fin/etat: debut et fin de la seq. dont l'etat est etat
// Match: resume des hits EST
// cons/incons: retour des valeurs
void CheckConsistency(int debut, int fin, int etat, 
		      unsigned char *Match, int * cons, int* incons);

// Verif coherence EST: calcul le nombre de nuc. coherents et
// incoherents avec les match est
// debut/fin/etat: debut et fin de la seq. dont l'etat est etat
// cons/incons: retour des valeurs
// WARNING : A modifier, utilise ESTMATCH_TMP (Cf struct DATA) !!!!!!!!!!!!
void CheckConsistency2(int debut, int fin, int etat, 
		       int * cons, int* incons, DNASeq *);

// Tdebut/fin = debut/fin de transcript
// debut/fin = debut/fin de traduit
void ESTSupport(char * Choice, int Tdebut, int Tfin, int debut,
		int fin,  Hits **HitTable, int Size);
