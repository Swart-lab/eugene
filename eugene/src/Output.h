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
#include "Prediction.h"

void Output (DNASeq *, Prediction *, int sequence, int argc, char * argv[]);

// Convertit les phases 0-6 en 1 2 3 -1 -2 -3 0
int PhaseAdapt(char);

// Verif coherence EST: calcul le nombre de nuc. coherents et
// incoherents avec les match est
// debut/fin/etat: debut et fin de la seq. dont l'etat est etat
// cons/incons: retour des valeurs
// WARNING : A modifier, utilise ESTMATCH_TMP (Cf struct DATA) !!!!!!!!!!!!
void CheckConsistency(int debut, int fin, int etat, 
		      int * cons, int* incons, DNASeq *);
