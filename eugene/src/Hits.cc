/*****************************************************************************/
/*             Copyright (c) 2004 by INRA. All rights reserved.              */
/*                 Redistribution is not permitted  without                  */
/*                 the express written permission of  INRA.                  */
/*                   Mail : eugene@ossau.toulouse.inra.fr                    */
/*---------------------------------------------------------------------------*/
/* File         : EuGeneTk/EuGene/Hits.cc                                    */
/* Description  : Definitions for a class representing alignements with      */
/*                genomic seq.                                               */
/* Authors      : P.Bardou, S.Foissac, M.J.Cros, A.Moisan, T.Schiex          */
/* History      : July 2004                                                  */
/*****************************************************************************/

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#include "Const.h"
#include "Hits.h"
#include <ctype.h>
#include <stdlib.h>
#include <assert.h>
#ifdef STDC_HEADERS
#include <string.h>
#else
#include <strings.h>
#endif
#include "System.h"

// ---------------------------------------------------------------------
//  Default constructor.
// ---------------------------------------------------------------------
Block :: Block ()
{
  Prev = Next = NULL;
}

// ---------------------------------------------------------------------
//  Constuct from Deb/Fin
// ---------------------------------------------------------------------
Block ::  Block(int start, int end,int lstart,int lend,int Ph,int Scr)
{
  Start  = start;
  End    = end;
  LStart = lstart;
  LEnd   = lend;
  Phase  = Ph;
  Score  = Scr;
  Prev   = Next = NULL;
}

// ---------------------------------------------------------------------
//  Insere apres
// ---------------------------------------------------------------------
void Block ::  AddBlockAfter(int start,int end,int lstart,int lend,int Ph,int Scr)
{
  Block *ABlock = new Block(start,end,lstart,lend,Ph,Scr);
  this->Next    = ABlock;
  ABlock->Prev  = this;
}

// ---------------------------------------------------------------------
//  Destroy this Block and all blocks after him
// ---------------------------------------------------------------------
Block :: ~ Block  ()
{
  delete Next;
}

// ---------------------------------------------------------------------
//  Default constructor.
// ---------------------------------------------------------------------
Hits :: Hits  ()
{
  Rejected = 0;
  Start    = 0;
  End      = 0;
  Evalue   = 0.0;
  Level    = 0;
  Support  = 0;
  Length   = 0;
  Strand   = 0;
  NGaps    = 0;
  Name     = NULL;
  Match    = NULL;
  Next     = NULL;
}

// ---------------------------------------------------------------------
//  Construct from a char*...
// ---------------------------------------------------------------------
Hits :: Hits  (char* name, int length, char strand, int deb, int fin, 
	       int ldeb, int lfin, int Ph, int Scr, double Prob, int level,
	       int sup)
{
  Rejected = 0;
  Name     = new char[strlen(name)+1];
  strcpy(Name,name);
  Strand = strand;
  Length = length;
  Start  = deb;
  End    = fin;
  Evalue = Prob;
  Level  = level;
  Support= sup;
  Match  = new Block(deb,fin,ldeb,lfin,Ph,Scr);
  NGaps  = 0;
  Next   = NULL;
}

// ---------------------------------------------------------------------
//  Read a table of Hits from a file
// ---------------------------------------------------------------------
Hits* Hits::ReadFromFile(FILE* HitFile, int *NumHits, int level, int margin) 
{
  char   *HitId, *PHitId;
  int    deb, fin, phase, Pphase, HSPDeb, HSPFin, poids, read;
  double evalue, Pevalue;
  char   A[128], B[128];
  Block *ThisBlock = NULL;
  Hits  *OneHit    = NULL, *ThisHit = this, *AllHit = this;
  const int MaxHitLen = 15000;

  Pevalue = -1.0;
  Pphase = 0;
  A[0]    = B[0] = 0;
  HitId   = A;
  PHitId  = B;
  
  if (ThisHit != NULL)
    for (int i=0; i<*NumHits-1; i++) ThisHit = ThisHit->Next;

  while ((read=fscanf(HitFile,"%d %d %d %lf %d %s %d %d\n", &deb, &fin, 
		      &poids, &evalue, &phase, HitId, &HSPDeb, &HSPFin)) == 8)
    {
      if (phase < 0) {
	int tmp = deb;
	deb     = fin;
	fin     = tmp;
	tmp     = HSPDeb;
	HSPDeb  = HSPFin;
	HSPFin  = tmp;
      }

      if (abs(fin-deb) > MaxHitLen) {
	fprintf(stderr,"Similarity of extreme length rejected. Check %s\n",
		HitId);
	continue;
      }

      if ((strcmp(HitId,PHitId) == 0)      && (phase*Pphase >= 0) &&
	  (deb + margin > ThisBlock->End)  &&
	  (phase >= 0  &&  (HSPDeb + margin > ThisBlock->LEnd)    ||
	   phase < 0   &&  (HSPDeb - margin < ThisBlock->LEnd)))
	// si HitId et PHitId sont égaux, alors il y a un Hit en cours
	// de meme nom on verifie que c'est bien compatible en terme
	// de position (sur l'est et le genomique) et en e-value, en 
	// prenant en compte le brin !
	{
	  ThisHit->NGaps++;
	  ThisHit->End = fin-1;
	  ThisBlock->AddBlockAfter(deb-1,fin-1,HSPDeb,HSPFin,phase,poids);
	  ThisBlock = ThisBlock->Next;
	}
      else {
	(*NumHits)++;
	OneHit = new Hits(HitId, poids, (phase>0 ? 1 : -1), deb-1, fin-1,
			  HSPDeb, HSPFin, phase, poids, evalue, level, 0);
	ThisBlock = OneHit->Match;
	PHitId    = HitId;
	Pevalue   = evalue;
	Pphase = phase;

	if (AllHit == NULL) {
	  AllHit = OneHit;
	  ThisHit = OneHit;
	}
	else {
	  ThisHit->Next = OneHit;
	  ThisHit       = OneHit;
	}
      }
    }
  if (read != EOF) fprintf(stderr,"Incorrect similarity file !\n");

  return AllHit;
}

// ---------------------------------------------------------------------
//  Destroy this alignement
// ---------------------------------------------------------------------
Hits :: ~ Hits  ()
{
  delete Match;
  delete [] Name;
  delete Next;
}
