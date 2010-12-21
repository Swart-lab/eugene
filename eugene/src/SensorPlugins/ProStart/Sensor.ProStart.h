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
// File:     Sensor.ProStart.h
// Contents: Sensor ProStart
// ------------------------------------------------------------------

#ifndef  SENSOR_PROSTART_H_INCLUDED
#define  SENSOR_PROSTART_H_INCLUDED

#include "../../Sensor.h"

#define REAL float
#define DEFAULTRBSENERGY 4.0
#define DEFAULT_POOR_ENERGY 4.0;

const double       MAXRBSPENALTY = -50.0;
const double       RT            = ( 1.942e-3 ) * ( 300.0 ); // R for kcal/mol, 300 K
const unsigned int MAXRBSLEN     = 12;
const int          MAXFREELEN    = 3;


// ************
// AlVertex aligment graph vertex
// ************
class AlVertex
{
	private:

	public:
		// Cost in case of
		//  - 0: match (helix)      Pair
		//  - 1: bulge (cote rRNA)  BulgeR
		//  - 2: bulge (cote mRNA)  BulgeM
		//  - 3: internal loop      IntLoop
		float E[4];
		// cas de bulge (cote rRNA), longueur du bulge
		short BulgeRLen;
		// cas de bulge (cote mRNA) et la longueur du bulge
		short BulgeMLen;

		// cas de boucle interne et les deux longueurs
		short IntLoopRLen, IntLoopMLen;

		// Memorisation du chemin
		// - 0: le noeud (0,1,2,3)
		// - 1: le delta en i
		// - 1: le delta en j
		short Parent[4][3];

		AlVertex();
		AlVertex (short, short, short, short, float, float, float, float);
		~AlVertex();
};


/*************************************************************
 **                     SensorProStart                      **
 *************************************************************/
class SensorProStart : public Sensor
{
	private:
		std::vector<double> vScoreF;       // vScoreF [i] gets the score at the position i on forward strand
		std::vector<double> vScoreR;       // vScoreR [i] gets the score at the position i on reverse strand
		std::vector<bool>   vDegeneratedF; // vDegeneratedF[i] = true means the codon start on forward, at the position i is degenerated
		std::vector<bool>   vDegeneratedR; // vDegeneratedR[i] = true means the codon start on reverse, at the position i is degenerated

		vector<vector<AlVertex> > AlignGraph;
		// ------------------------------------------
		// Stacking energy. Index[i][j][k][l]  is for 5' >-IJ-> 3'
		// ------------------------------------------ 3' <-KL<- 5'
		float  StackEnergy[4][4][4][4];
		float  IntLoopEnergy[30];
		float  BulgeEnergy[30];
		char   RBS[MAXRBSLEN];
		double alpha;
		double beta;
		int    matchlen;
		int    matchoffset;

		void SearchProStart (DNASeq *X);
		void ReadScore();
		REAL  StackingEnergy (char F1, char F2, char R1, char R2);
		int  Nuc2Index       (char nuc);
		REAL LoopPenalty     (int L, int R);
		REAL BulgePenalty    (int L);
		int  Pair            (char P, char Q );
		REAL Fold            (char * R, char * M );

	public:
		SensorProStart ( int n, DNASeq *X );
		virtual ~SensorProStart ();
		virtual void Init (DNASeq * );
		virtual void GiveInfo (DNASeq *X, int, DATA * );
		virtual void Plot (DNASeq *X );
		virtual void PostAnalyse (Prediction *, FILE * );
};

extern "C" SensorProStart * builder0 ( int n, DNASeq *X ) {  return new SensorProStart ( n, X ); }


#endif
