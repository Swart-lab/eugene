// ------------------------------------------------------------------
// Copyright (C) 2014 INRA <eugene-help@mulcyber.toulouse.inra.fr>
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
// File:     WAMbuilder.cc
//
// This program build a Weight Array Model (WAM) for a WAM eugene sensor.
// The goal is to model a signal in a genomic sequence (splice site, start)
// Build the model from a multifasta file containing sequences around the
// around the signal. Save results in binary files.
//
// USAGE: egn_WAMBuilder motifFile pat uc dc markovOrder outfilePrefix
//  motifFile     multifasta file of the motifs,
//  pat           signal pattern,
//  uc            upstream sequence length,
//  dc            downstream length, 
//  markovOrder   markovian order (start at 0),
//  outfilePrefix prefix of the results binary files.
// EXAMPLES: 
// egn_WAMBuilder ARA.don.TP.fa GT 4 4 0 WAM.don.order0.TP.
// egn_WAMBuilder ARA.don.TP.fa GT 4 4 1 WAM.don.order1.TP.
// egn_WAMBuilder ARA.acc.TN.fa AG 4 4 0 WAM.acc.order0.TN.
//
// ------------------------------------------------------------------

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <vector>

#include "../markov.h"
#include "../markov.cc"

// -------------------------------------------------------------------------
//            MAIN
// -------------------------------------------------------------------------
int main(int argc, char *argv[])
{
    if (argc != 7)
    {
        fprintf(stderr, "\
This program build a Weight Array Model (WAM) for a WAM eugene sensor. \n\
The goal is to model a signal in a genomic sequence (splice site, start) \n\
Build the model from a multifasta file containing sequences around the \n\
around the signal. Save results in binary files. \n\
\nUSAGE: %s motifFile pat uc dc markovOrder outfilePrefix \n\
  motifFile     multifasta file of the motifs, \n\
  pat           signal pattern, (example: GT)\n\
  uc            upstream sequence length, \n\
  dc            downstream length, \n\
  markovOrder   markovian order (start at 0), \n\
  outfilePrefix prefix of the result binary files.\
\n\nEXAMPLES: \
\n%s ARA.don.TP.fa GT 4 4 0 WAM.don.order0.TP. \
\n%s ARA.don.TP.fa GT 4 4 1 WAM.don.order1.TP. \
\n%s ARA.acc.TN.fa AG 4 4 0 WAM.acc.order0.TN.\n", argv[0], argv[0], argv[0], argv[0]);
        exit(1);
    }

    int DEBUG = 0;
	int   i;
	FILE* file;
	
    char* motifFile  = argv[1];       // mutfifasta file
    char* signal     = argv[2];       // GC or GC or AT or ATG, etc
    int uc           = atoi(argv[3]); // upstream region length
    int dc           = atoi(argv[4]); // downstream region length
    int markovOrder  = atoi(argv[5]); // 0, 1, etc
    char* outname    = argv[6];


    int signalLen   = strlen(signal);
    int motifLen    = uc + signalLen + dc;
    int WAMlen      = motifLen - markovOrder;
    int basenamelen = strlen(outname);
    ChaineADN* ADN  = new ChaineADN;

	// Contain the word occurences in the input sequences
    std::vector<TabChaine<ChaineADN, int>*> WAMCOUNT; 
	// Contain probabilities 
    std::vector<TabChaine<ChaineADN, unsigned short int>*> WAMMOD;

	// vectors initialisation
    for (i = 0; i < WAMlen; i++)
    {
        WAMCOUNT.push_back(new TabChaine<ChaineADN, int>(markovOrder, ADN));
        WAMMOD.push_back(new TabChaine<ChaineADN, unsigned short int>(markovOrder, ADN));
    }

    //----------------------------------------------------------------//
    // COUNT words occurences for each position
    fprintf(stderr, " - counting words occurences (order %d)  ... position...", markovOrder);

    for (i = 0; i < WAMlen; i++)
    {
        file = fopen(motifFile, "rt");

        if (file == NULL)
        {
            printf("Cannot open consensus sequence file %s\n", motifFile);
            exit(1);
        }

        fprintf(stderr, "%d...", i);
        WAMCOUNT[i]->fichier2compte(file, markovOrder + i, markovOrder + i);
        fclose(file);
    }

    fprintf(stderr, "done\n");

    if (DEBUG)
    {
        for (i = 0; i < WAMlen; i++)
        {
            printf("Position %d:\n", i);
            WAMCOUNT[i]->affichage(0);
        }
    }


    //----------------------------------------------------------------//
    // ADD pseudocount to each word occurence (to avoid null probs)
    fprintf(stderr, " - adding pseudocounts...");

    for (i = 0; i < WAMlen; i++)
    {
        WAMCOUNT[i]->pseudocount();
    }

    fprintf(stderr, "done\n");

    if (DEBUG)
    {
        for (i = 0; i < WAMlen; i++)
        {
            printf("Position %d:\n", i);
            WAMCOUNT[i]->affichage(0);
        }
    }


    //----------------------------------------------------------------//
    // COMPUTE PROBS (store results in matrix WAMMOD)
    fprintf(stderr, " - computing probs from frequences for each position...");

    for (i = 0; i < WAMlen; i++)
    {
        WAMMOD[i]->compte2probas(WAMCOUNT[i]);
        fprintf(stderr, "%d...", i);
    }

    fprintf(stderr, "done\n");

    if (DEBUG)
    {
        for (i = 0; i < WAMlen; i++)
        {
            printf("position %d:\n", i);
            WAMMOD[i]->affichage(0);
        }
    }

    //----------------------------------------------------------------//
    // SAVE data in the matrix file
    fprintf(stderr, " - saving data in matrix files... %s...", outname);

    for (i = 0; i < WAMlen; i++)
    {
        char* name = new char[basenamelen + 4];
        name[basenamelen + 3] = '\0';
        sprintf(name, "%s", outname);
		sprintf(name + basenamelen, "%02d", i);

        file = fopen(name, "wb");

        if (file == NULL)
        {
            printf("Cannot open matrix file %s\n", name);
            exit(1);
        }

        WAMMOD[i]->sauve2fichier(file);
        fclose(file);
        fprintf(stderr, "%d...", i);
        delete [] name;
    }

    fprintf(stderr, "done\n");

    for (i = 0; i < WAMlen; i++)
    {
        WAMCOUNT[i]->~TabChaine();
        WAMMOD[i]->~TabChaine();
    }

    delete ADN;

}
