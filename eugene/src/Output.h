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
// File:     Output.h
// Contents: functions to make the outputs
// ------------------------------------------------------------------

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <sstream>
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

void Output (DNASeq *, MasterSensor *, Prediction *, int sequence, 
	     int argc, char * argv[], FILE* f);

// -ph print on stdout the begin of the HTML output
void StartHTML(char* html_dir);

// -ph print on stdout the end of the HTML output
void EndHTML();

// Convertion int en string
std::string to_string(int);

// Convertit les phases 0-6 en 1 2 3 -1 -2 -3 0
int PhaseAdapt(char);

// Verif coherence EST: calcul le nombre de nuc. coherents et
// incoherents avec les match est
// debut/fin/etat: debut et fin de la seq. dont l'etat est etat
// cons/incons: retour des valeurs
void CheckConsistency(int debut, int fin, int etat, 
		      int * cons, int* incons, DNASeq *, MasterSensor *);
