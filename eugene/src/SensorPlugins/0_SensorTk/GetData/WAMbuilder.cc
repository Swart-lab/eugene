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
// Usage: %s motifFile species site pat seqType uc dc markovOrder [wamDir]
// motifFile                            Multifasta file of the motifs,
// species                              Species
// site          [donor|acceptor|start] Signal type
// pat                                  Signal pattern
// seqType       [TP|TN]                True positive or true negative sequences
// uc                                   Upstream sequence length
// dc                                   Downstream sequence length
// markovOrder                          Markovian order (starts at 0)
// 
// Options:
// wamDir                               Directory where save the wam models. (Default $EUGENEDIR/models/WAM)
// EXAMPLES: 
// egn_WAMbuilder ARA.don.TP.fa arath donor    GT TP 60 60 0
// egn_WAMbuilder ARA.don.TP.fa arath donor    GT TP 60 60 1
// egn_WAMbuilder ARA.acc.TN.fa arath acceptor AG TN 60 60 0 

//
// ------------------------------------------------------------------

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <vector>
#include <dirent.h>

#include <sys/types.h>
#include <sys/stat.h>
#include "../markov.h"
#include "../markov.cc"



int CreateDirectory(const char* path)
{
    DIR *dir;
    dir = opendir(path);
    if (!dir)
    {
        int statusMkdir = mkdir (path, 0711);
        if (statusMkdir != 0)
        {
            printf ("ERROR: cant create the directory %s.\n", path);
            return (-1); 
        }
        else
        {
            printf ("Directory %s created.\n", path);
        }  
    }
    else
    {
        printf ("Directory %s already exists.\n", path);
        closedir(dir);
    }
    
    return 0;
}

void PrintHelp (char* pg)
{
    fprintf(stderr, "\
    This program builds a Weight Array Model (WAM) for a WAM eugene sensor. \n\
    The goal is to model a signal in a genomic sequence (splice site, start).\n\
    Build the model from a multifasta file containing sequences around the \n\
    signal. Save results in binary files. \n\
    \nUsage: %s motifFile species site pat seqType uc dc markovOrder [wamDir] \n\
    motifFile                            Multifasta file of the motifs, \n\
    species                              Species\n\
    site          [donor|acceptor|start] Signal type\n\
    pat                                  Signal pattern\n\
    seqType       [TP|TN]                True positive or true negative sequences\n\
    uc                                   Upstream sequence length\n\
    dc                                   Downstream sequence length\n\
    markovOrder                          Markovian order (starts at 0)\n\
    \nOptions: \n\
    wamDir                               Directory where save the wam models. (Default $EUGENEDIR/models/WAM)\n\
    \n\nEXAMPLES: \
    \n%s ARA.don.TP.fa arath donor    GT TP 60 60 0\
    \n%s ARA.don.TP.fa arath donor    GT TP 60 60 1\
    \n%s ARA.acc.TN.fa arath acceptor AG TN 60 60 0\n", pg, pg, pg, pg);
    
    exit(1); 
}


// -------------------------------------------------------------------------
//            MAIN
// -------------------------------------------------------------------------
int main(int argc, char *argv[])
{
    
    int DEBUG = 0;
    int   i;
    FILE* file;
    
    
    if (argc < 9)
    {
        PrintHelp(argv[0]);
    }
    
    char* outRoot = new char[FILENAME_MAX+1];
    if (argc == 9)
    {
        char* eugenedir = getenv("EUGENEDIR");
        if (eugenedir == NULL) 
        {
            fprintf(stderr, ">ERROR: the parameter 'wamDir' should be defined<\n");
            PrintHelp(argv[0]);
        }
        sprintf(outRoot, "%s/models/WAM", eugenedir);
    }
    else
    {
        sprintf(outRoot, argv[9]);
    }

	
    char* motifFile  = argv[1];       // multifasta file
    char* species    = argv[2];       // species
    char* site       = argv[3];       // acceptor donor or start
    char* signal     = argv[4];       // GC or GC or AT or ATG, etc
    char* seqType    = argv[5];       // TP or TN
    int uc           = atoi(argv[6]); // upstream region length
    int dc           = atoi(argv[7]); // downstream region length
    int markovOrder  = atoi(argv[8]); // 0, 1, etc
   
    
    if ( strcmp(seqType, "TP") != 0 && strcmp(seqType, "TN") != 0 )
    {
        fprintf (stderr, ">ERROR: allowed values for seqType are TP and TN<\n\n");
        PrintHelp(argv[0]);  
    }
    
    if ( strcmp(site, "donor") != 0 && strcmp(site, "acceptor") != 0 && strcmp(site, "start") != 0 )
    {
        fprintf (stderr, ">ERROR: allowed values for site are: donor, acceptor and start<\n\n");
        PrintHelp(argv[0]);
  
    }
    if (markovOrder < 0)
    {    
        fprintf (stderr, ">ERROR: markovOrder has to be equal or greater than 0<\n\n");
        PrintHelp(argv[0]);
            
    }
    assert(uc > 0);
    assert(dc > 0);
    
    // Check that output directory exists
    DIR *rootDir = opendir(outRoot);
    if (!rootDir)
    {
        fprintf (stderr, ">ERROR: %s is not an existing directory<\n", outRoot);
        exit(1);
    }
    closedir(rootDir);
    
    // Create the species directory
    char* speciesPath = new char[FILENAME_MAX+1];
    sprintf(speciesPath, "%s/%s/", outRoot, species);
    int statusMkdir = CreateDirectory(speciesPath);
    if (statusMkdir != 0)
    {
        exit(1); 
    }
    
    // Create the site directory
    char* sitePath = new char[FILENAME_MAX+1];
    sprintf(sitePath, "%s/%s/", speciesPath, site);
    statusMkdir = CreateDirectory(sitePath);
    if (statusMkdir != 0)
    {
        exit(1); 
    }
    // build outfile root name
    // Example: EUGENEDIR/models/WAM/plants/acceptor/WAM.o04.uc60.dc60.AG.TP.
    char *outname = new char[FILENAME_MAX+1];
    sprintf(outname, "%s/WAM.o%02d.uc%03d.dc%03d.%s.%s.", sitePath, markovOrder, uc, dc, signal, seqType);
    
    
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

    // check the length od the motif sequences
    file = fopen(motifFile, "rt");

    if (file == NULL)
    {
        printf("Cannot open consensus sequence file %s\n", motifFile);
        exit(1);
    }
    char* Sequence=NULL;

    while (WAMCOUNT[0]->fichier2seq(file, Sequence))
    {
        int seqLen = strlen(Sequence);

        if (seqLen != motifLen)
        {
            printf("Error reading >%s<: sequence \"%s\" length is different to expected length (%d nt)\n", motifFile, Sequence, motifLen);
            exit(1);
        }

        // check that the sequence includes the signal
        std::string SequenceString = Sequence;
        if (SequenceString.compare(uc, signalLen, signal) != 0)
        {
            printf("Error reading >%s<: the sequence \"%s\" does not contain \"%s\" at the position %d.\n", motifFile, Sequence, signal, (uc+1));
            exit(1);
        }
        
        free(Sequence);
        Sequence = NULL;
    };

    fclose(file);
	
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
    
    int digitnb = (int) (1 + log10(WAMlen-1)); // nombre de chiffres de l'entier (WAMlen-1)
    
    for (i = 0; i < WAMlen; i++)
    {
        char* name = new char[basenamelen + digitnb + 2];
        name[basenamelen + digitnb +1 ] = '\0';
        sprintf(name, "%s", outname);
		sprintf(name + basenamelen, "%0*d", digitnb, i);

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


