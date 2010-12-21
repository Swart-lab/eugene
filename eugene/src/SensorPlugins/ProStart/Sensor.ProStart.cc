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
// File:     Sensor.ProStart.cc
// Contents: Sensor ProStart
// ------------------------------------------------------------------

#include "Sensor.ProStart.h"

extern Parameters PAR;

#define NORM(x,n) Max(3.0,(x)+(n))/(n)

/*************************************************************
 **                      AlVertex                           **
 *************************************************************/
// ----------------------
//  Default constructor.
// ----------------------
AlVertex::AlVertex()
{
}

// ----------------------
//  Constructor.
// ----------------------
AlVertex::AlVertex ( short brl, short bml, short ilrl, short ilml, float e1, float e2, float e3, float e4 )
{
	BulgeRLen   = brl;
	BulgeMLen   = bml;
	IntLoopRLen = ilrl;
	IntLoopMLen = ilml;
	E[0] = e1;
	E[1] = e2;
	E[2] = e3;
	E[3] = e4;
}

// ----------------------
//  Destructor
// ----------------------
AlVertex::~AlVertex()
{
}


/*************************************************************
 **                      SensorProStart                       **
 *************************************************************/
// ----------------------
//  Default constructor.
// ----------------------
SensorProStart :: SensorProStart ( int n, DNASeq *X ) : Sensor ( n )
{
	vScoreF.clear();
	vScoreR.clear();
	vDegeneratedF.clear();
	vDegeneratedR.clear();

	// Get the parameters from the parameter file
	matchlen      = PAR.getD ( "ProStart.matchlen" );
	matchoffset   = PAR.getD ( "ProStart.matchoffset" );
	strcpy ( RBS, PAR.getC("ProStart.RBSPattern"));

	AlignGraph.resize ( 1+MAXRBSLEN, vector<AlVertex> ( 1 + matchlen ) );

	ReadScore(); // Read the energy scores
}


// ----------------------
//  Default destructor.
// ----------------------
SensorProStart :: ~SensorProStart ()
{
	vScoreF.clear();
	vScoreR.clear();
	vDegeneratedF.clear();
	vDegeneratedR.clear();
}

// ----------------------
//  Init .
// ----------------------
void SensorProStart :: Init ( DNASeq *X )
{
	vScoreF.clear();
	vScoreR.clear();
	vDegeneratedF.clear();
	vDegeneratedR.clear();

	// Get the optimizable parameters
	alpha = PAR.getD ( "ProStart.alpha*" );
	beta  = exp ( PAR.getD ( "ProStart.beta*" ) );

	SearchProStart(X); // Search the start codons in the sequence

	if ( PAR.getI ( "Output.graph" ) ) Plot (X);
}


// ------------------------------------------------------------------------
// Read Energy scores from Zucker/Turner.  The stack.dat file must be
// modified as follows: all non numerical lines are removed and "."
// are replaced by 1e10 or 0.8 (min energy for small loops).
// ------------------------------------------------------------------------
void SensorProStart::ReadScore ()
{
	FILE *Data;
	double e, f;
	int    i, j, k, l;
	char *models_dir = new char[FILENAME_MAX+1];

	strcpy ( models_dir, PAR.getC ( "eugene_dir" ) );
	strcat ( models_dir, MODELS_DIR );
	if ( ! ( Data = FileOpen ( models_dir, PAR.getC ( "ProStart.stackFile" ), ( char * ) "r" ) ) )
	{
		fprintf ( stderr, "cannot open RNA binding energy parameter file ProStart.stackFile\n" );
		exit ( 2 );
	}

	for ( i=0; i<4; i++ )
		for ( j=0; j<4; j++ )
			for ( k=0; k<4; k++ )
				for ( l=0; l<4; l++ )
				{
					fscanf ( Data,"%lf",&e );
					StackEnergy[i][j][k][l] = ( ( e >= 9e9 ) ? INFINITY : e );
				}
	fclose ( Data );

	if ( ! ( Data = FileOpen ( models_dir, PAR.getC ( "ProStart.loopFile" ), ( char * ) "r" ) ) )
	{
		fprintf ( stderr, "cannot open RNA binding energy parameter file ProStart.loopFile\n" );
		exit ( 2 );
	}

	for ( i = 0; i < 30; i++ )
	{
		fscanf ( Data, "%*d %lf %lf %*f\n", &e, &f );
		IntLoopEnergy[i] = e;
		BulgeEnergy[i]   = f;
	}

	fclose ( Data );
}


// ------------------------------------------------------------------------
// Convert nuc to 1/4 according to zucker/Turner order (ACGU)
// ------------------------------------------------------------------------
int SensorProStart::Nuc2Index ( char nuc )
{
	switch ( nuc )
	{
		case 'a': return 0;
		case 'c': return 1;
		case 'g': return 2;
		case 't': return 3;
		default : fprintf ( stderr, "ERROR:  Bad nucleotide %c\n", nuc );
	}
}

// ------------------------------------------------------------------------
// Energy penalty for an internal loop with L/R bases on each side
// ------------------------------------------------------------------------
REAL SensorProStart::LoopPenalty ( int L, int R )
{
	int  TheMin;

	if ( L <= 0 || R <= 0 )
		return  INFINITY;

	TheMin = Min ( L,R );

	if ( TheMin > 30 ) TheMin = 30;

	//apply Ninio rule: 0.5 penalty per assymetry with a maximum of 3.0
	return IntLoopEnergy[TheMin-1] + Min ( 3.0,abs ( L-R ) *0.5 );
}

// ------------------------------------------------------------------------
// Energy penalty for a bulge with L bases
// ------------------------------------------------------------------------
REAL  SensorProStart::BulgePenalty ( int L )
{
	if ( L <= 0 )
		return  INFINITY;

	if ( L > 30 ) L=30;

	return  BulgeEnergy[L-1];
}


// ------------------------------------------------------------------------
// Possible match (Watson/Crick or Wobble)
// ------------------------------------------------------------------------
int SensorProStart::Pair ( char P, char Q )
{
	switch ( P )
	{
		case  'a' :
			return ( Q == 't' );
		case  'c' :
			return ( Q == 'g' );
		case  'g' :
			return ( Q == 'c' || Q == 't' );
		case  't' :
			return ( Q == 'a' || Q == 'g' );
	}
	return  0;
}

// ------------------------------------------------------------------------
// Stacking energy according to Zucker/Turner web site at 37 C in      F1F2
// kcal/mol.F1/2 (5' -> 3', 1st strand) is bound with R1/2 (3' -> 5')  R1R2
// ------------------------------------------------------------------------
inline REAL SensorProStart::StackingEnergy ( char F1, char F2, char R1, char R2 )
{
	return StackEnergy[Nuc2Index ( F1 ) ][Nuc2Index ( F2 ) ][Nuc2Index ( R1 ) ][Nuc2Index ( R2 ) ];
}

// ------------------------------------------------------------------------
// Returns the optimal RNA free energy interaction of string R (rRNA)
// with string M (mRNA). Ignore the constant energy difference for 2
// molecules
// ------------------------------------------------------------------------
REAL SensorProStart::Fold ( char * R, char * M )
{
	REAL  Best, TheBest, Optimum, Value;
	int  a, b, i, j, LR, LM, LenR, LenM;

	LR = strlen ( R );
	LM = strlen ( M );

	for ( i = 0;  i < LR;  i ++ )
		R[i] = tolower ( R[i] );
	for ( j = 0;  j < LM;  j ++ )
		M[i] = tolower ( M[i] );
	for ( i = 1;  i <= LR;  i ++ )
		AlignGraph[i][0] = AlVertex ( i, 0, 0, 0, INFINITY, 0.0, INFINITY, INFINITY );
	for ( j = 1;  j <= LM;  j ++ )
		AlignGraph[0][j] = AlVertex ( 0, j, 0, 0, INFINITY, INFINITY, 0.0, INFINITY );
	AlignGraph[0][0] = AlVertex ( 0, 0, 0, 0, INFINITY, 0.0, 0.0, 0.0 );


	Optimum = INFINITY;
	for ( j = 1;  j <= LM;  j ++ )
	{
		TheBest = INFINITY;

		for ( i = 1;  i <= LR;  i ++ )
		{
			// -----------------------------------------------
			// Les bulges cote rRNA
			// -----------------------------------------------
			Best = INFINITY;
			LenR = 0;

			// A partir d'une paire, on ouvre un bulge de long. 1
			Value = AlignGraph[i-1][j].E[0];
			if ( Value < Best )
			{
				Best = Value;
				LenR = 1;
				AlignGraph[i][j].Parent[1][0] = 0;
				AlignGraph[i][j].Parent[1][1] = -1;
				AlignGraph[i][j].Parent[1][2] = 0;
			}

			// A partir d'un bulge rRNA, on l'etend
			Value = AlignGraph[i-1][j].E[1];
			a = AlignGraph[i-1][j].BulgeRLen;
			if ( Value < Best && a < MAXFREELEN )
			{
				Best = Value;
				LenR = 1+a;
				AlignGraph[i][j].Parent[1][0] = 1;
				AlignGraph[i][j].Parent[1][1] = -1;
				AlignGraph[i][j].Parent[1][2] = 0;
			}

			// On prend le meilleur des deux
			AlignGraph[i][j].E[1]      = Best;
			AlignGraph[i][j].BulgeRLen = LenR;

			TheBest = Min ( Best,TheBest );

			// -----------------------------------------------
			// Les bulges cote mRNA
			// -----------------------------------------------

			Best = INFINITY;
			LenM = 0;

			// A partir d'une paire, ca ouvre un bulge de long. 1
			Value = AlignGraph[i][j-1].E[0];
			if ( Value < Best )
			{
				Best = Value;
				LenM = 1;
				AlignGraph[i][j].Parent[2][0] = 0;
				AlignGraph[i][j].Parent[2][1] = 0;
				AlignGraph[i][j].Parent[2][2] = -1;
			}

			// A partir d'un bulge, ca l'etend
			Value = AlignGraph[i][j-1].E[2];
			a = AlignGraph[i][j-1].BulgeMLen;
			if ( Value < Best && a < MAXFREELEN )
			{
				Best = Value;
				LenM = a+1;
				AlignGraph[i][j].Parent[2][0] = 2;
				AlignGraph[i][j].Parent[2][1] = 0;
				AlignGraph[i][j].Parent[2][2] = -1;
			}

			// On prend le meilleur des deux
			AlignGraph[i][j].E[2]      = Best;
			AlignGraph[i][j].BulgeMLen = LenM;

			// -----------------------------------------------
			// Les internal loops
			// -----------------------------------------------

			Best = INFINITY;
			LenR = LenM = 0;

			// A partir d'une paire, on ouvre une internal loop (1/1)
			Value = AlignGraph[i-1][j-1].E[0];
			if ( Value < Best )
			{
				Best = Value;
				LenM = LenR = 1;
				AlignGraph[i][j].Parent[3][0] = 0;
				AlignGraph[i][j].Parent[3][1] = -1;
				AlignGraph[i][j].Parent[3][2] = -1;
			}

			// A partir d'un bulge mRNA, on fait une internal loop (1/L)
			Value = AlignGraph[i-1][j].E[2];
			if ( Value < Best )
			{
				Best = Value;
				LenR = 1;
				LenM = AlignGraph[i-1][j].BulgeMLen;
				AlignGraph[i][j].Parent[3][0] = 2;
				AlignGraph[i][j].Parent[3][1] = -1;
				AlignGraph[i][j].Parent[3][2] = 0;
			}

			// A partir d'un bulge rRNA, on fait une internal loop (L/1)
			Value = AlignGraph[i][j-1].E[1];
			if ( Value < Best )
			{
				Best = Value;
				LenR = AlignGraph[i][j-1].BulgeRLen;
				LenM = 1;
				AlignGraph[i][j].Parent[3][0] = 1;
				AlignGraph[i][j].Parent[3][1] = 0;
				AlignGraph[i][j].Parent[3][2] = -1;
			}

			// A partir d'une internal loop, on l'allonge cote rRNA seulement
			Value = AlignGraph[i-1][j].E[3];
			a = AlignGraph[i-1][j].IntLoopRLen;
			if ( Value < Best && a < MAXFREELEN )
			{
				Best = Value;
				LenR = a+1;
				LenM = AlignGraph[i-1][j].IntLoopMLen;
				AlignGraph[i][j].Parent[3][0] = 3;
				AlignGraph[i][j].Parent[3][1] = -1;
				AlignGraph[i][j].Parent[3][2] = 0;
			}

			// A partir d'une internal loop, on l'allonge cote mRNA seulement
			Value = AlignGraph[i][j-1].E[3];
			a = AlignGraph[i][j-1].IntLoopMLen;
			if ( Value < Best && a < MAXFREELEN )
			{
				Best = Value;
				LenR = AlignGraph[i][j-1].IntLoopRLen;
				LenM = a+1;
				AlignGraph[i][j].Parent[3][0] = 3;
				AlignGraph[i][j].Parent[3][1] = 0;
				AlignGraph[i][j].Parent[3][2] = -1;
			}

			// A partir d'une internal loop, on l'allonge des deux cotes
			Value = AlignGraph[i-1][j-1].E[3];
			a = AlignGraph[i-1][j-1].IntLoopRLen;
			b = AlignGraph[i-1][j-1].IntLoopMLen;
			if ( Value < Best && a < MAXFREELEN && b < MAXFREELEN )
			{
				Best = Value;
				LenR = a+1;
				LenM = b+1;
				AlignGraph[i][j].Parent[3][0] = 3;
				AlignGraph[i][j].Parent[3][1] = -1;
				AlignGraph[i][j].Parent[3][2] = -1;
			}

			// On prend le meilleur des differents cas
			AlignGraph[i][j].E[3] = Best;
			AlignGraph[i][j].IntLoopRLen = LenR;
			AlignGraph[i][j].IntLoopMLen = LenM;

			// -----------------------------------------------
			// Les helices
			// -----------------------------------------------
			Best = INFINITY;

			// si on ne peut pas lier, fini
			if ( !Pair ( R[i-1], M[j-1] ) )
				AlignGraph[i][j].E[0] = INFINITY;
			else
			{
				// on ferme un bulge cote ARN rib
				Value = AlignGraph[i-1][j-1].E[1]+
				        BulgePenalty ( AlignGraph[i-1][j-1].BulgeRLen );
				if ( Value < Best )
				{
					Best = Value;
					AlignGraph[i][j].Parent[0][0] = 1;
					AlignGraph[i][j].Parent[0][1] = -1;
					AlignGraph[i][j].Parent[0][2] = -1;
				}

				// on ferme un bulge cote ARN m
				Value = AlignGraph[i-1][j-1].E[2]+
				        BulgePenalty ( AlignGraph[i-1][j-1].BulgeMLen );
				if ( Value < Best )
				{
					Best = Value;
					AlignGraph[i][j].Parent[0][0] = 2;
					AlignGraph[i][j].Parent[0][1] = -1;
					AlignGraph[i][j].Parent[0][2] = -1;
				}

				// on ferme une boucle interne
				Value = AlignGraph[i-1][j-1].E[3]+
				        LoopPenalty ( AlignGraph[i-1][j-1].IntLoopRLen,
				                      AlignGraph[i-1][j-1].IntLoopMLen );
				if ( Value < Best )
				{
					Best = Value;
					AlignGraph[i][j].Parent[0][0] = 3;
					AlignGraph[i][j].Parent[0][1] = -1;
					AlignGraph[i][j].Parent[0][2] = -1;
				}

				// on prolonge une helice. Attention au debut des sequences
				// (R/M[i-2] indefinis => Value = INFINITY
				Value = AlignGraph[i-1][j-1].E[0];
				if ( Value != INFINITY ) Value += StackingEnergy ( R[i-2], R[i-1], M[j-2], M[j-1] );
				if ( Value < Best )
				{
					Best = Value;
					AlignGraph[i][j].Parent[0][0] = 0;
					AlignGraph[i][j].Parent[0][1] = -1;
					AlignGraph[i][j].Parent[0][2] = -1;
				}

				AlignGraph[i][j].E[0] = Best;
				TheBest = Min ( Best,TheBest );
			}
		}
		Optimum = Min ( TheBest,Optimum );
	}

	if ( Optimum > 100 ) Optimum = DEFAULT_POOR_ENERGY;

	return  Optimum;
}


// ------------------------
//  GiveInfo signal start.
// ------------------------
void SensorProStart::GiveInfo ( DNASeq *X, int pos, DATA *d )
{
	double score; // score de ProStart

	if (vScoreF[pos] != 0)
	{
		score = vScoreF[pos];
		// Apply the score to the signal
		d->sig[DATA::Start].weight[Signal::Forward] += log ( score);
		if 	(!vDegeneratedF[pos]) // not a degenerated start
			d->sig[DATA::Start].weight[Signal::ForwardNo] += log ( 1.0 - score );
	}

	if (vScoreR[pos] != 0)
	{
		score = vScoreR[pos];
		// Apply the score to the signal
		d->sig[DATA::Start].weight[Signal::Reverse] += log ( score );
		if 	(!vDegeneratedR[pos]) // not a degenerated start
			d->sig[DATA::Start].weight[Signal::ReverseNo] += log ( 1.0 - score );
	}
}


// ------------------------
//  Search start codon at each position of the sequence
// ------------------------
void SensorProStart :: SearchProStart (DNASeq *X)
{
	char   Buffer[1+matchlen];
	double StartDetected;
	double score;
	bool   isDegenerated;

	for (int pos = 0;  pos <= X->SeqLen;  pos++) 
	{
		score         = 0;
		isDegenerated = true;

		// strand +
		if ( (pos < X->SeqLen-2 ) && ( StartDetected = X->IsProStart ( pos, 1 ) ) != 0.0 )
		{
			if ( pos >= matchlen + matchoffset )
			{
				X->Transfer ( pos - matchlen - matchoffset, matchlen, Buffer, 1 );
				score = StartDetected * Max ( exp ( MAXRBSPENALTY ), X->StartPenalty ( pos, 1 ) *
											alpha/ ( 1.0 + beta * exp ( Fold ( RBS, Buffer ) /RT ) ) );
			}
			else
			{
				score = StartDetected * Max ( exp ( MAXRBSPENALTY ), X->StartPenalty ( pos, 1 ) *
											alpha/ ( 1.0 + beta*exp ( DEFAULTRBSENERGY/RT ) ) );
			}
			if ( StartDetected == 1 ) // not a degenerated start
				isDegenerated = false;

		}
		// Save the score and if the codon is degenerated
		vScoreF.push_back(score);
		vDegeneratedF.push_back(isDegenerated);

		// Reinit variables
		score         = 0;
		isDegenerated = true;
	
		// strand -
		if ( (pos > 2) && ( StartDetected = X->IsProStart ( pos-1,-1 ) ) != 0.0 )
		{
			if ( pos -1 + matchlen + matchoffset < X->SeqLen )
			{
				X->Transfer ( pos + matchlen + matchoffset -1, -matchlen, Buffer, 1 );
				score = StartDetected * Max ( exp ( MAXRBSPENALTY ), X->StartPenalty ( pos-1,-1 ) *
											alpha/ ( 1.0 + beta*exp ( Fold ( RBS, Buffer ) /RT ) ) );
			}
			else
			{
				score = StartDetected * Max ( exp ( MAXRBSPENALTY ), X->StartPenalty ( pos-1,-1 ) *
											alpha/ ( 1.0+beta*exp ( DEFAULTRBSENERGY/RT ) ) );
			}
			if ( StartDetected == 1 ) // not a degenerated start
				isDegenerated = false;
		}
		// Save the score and if the codon is degenerated
		vScoreR.push_back(score);
		vDegeneratedR.push_back(isDegenerated);
	}

	assert ( vScoreF.size()       == (X->SeqLen+1));
	assert ( vScoreR.size()       == (X->SeqLen+1));
	assert ( vDegeneratedF.size() == (X->SeqLen+1));
	assert ( vDegeneratedR.size() == (X->SeqLen+1));
}

// ----------------------------
//  Plot Sensor information
// ----------------------------
void SensorProStart :: Plot ( DNASeq *X )
{
	for (int pos=0; pos <= X->SeqLen; pos++) 
	{
		if (vScoreF[pos] != 0)
			PlotStart(pos, (pos%3)+1, NORM(log(vScoreF[pos]),10.0));

		if (vScoreR[pos] != 0)
			 PlotStart(pos, -((X->SeqLen-pos)%3)-1, NORM(log(vScoreR[pos]),10.0));
	}
}

// ------------------
//  Post analyse
// ------------------
void SensorProStart :: PostAnalyse ( Prediction *pred, FILE *MINFO )
{
}
