//
//   T. Schiex
//
//     File:  Hits.cc
//     Version:  1.0
//
//    Copyright (c) 2000 by Thomas Schiex All rights reserved.
//    Redistribution is not permitted without the express written
//    permission of the authors.
//
//  Definitions for a class representing alignements wrt to a seq.
//

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
  Prev =  Next = NULL;
}
// ---------------------------------------------------------------------
//  Constuct from Deb/Fin
// ---------------------------------------------------------------------
Block ::  Block(int start, int end,int lstart,int lend,int Ph,int Scr)
{
  Start = start;
  End = end;
  LStart = lstart;
  LEnd = lend;
  Phase = Ph;
  Score = Scr;
  Prev = Next = NULL;
}

// ---------------------------------------------------------------------
//  Insere apres
// ---------------------------------------------------------------------
void Block ::  AddBlockAfter(int start,int end,int lstart,int lend,int Ph,int Scr)
{
  Block * ABlock = new Block(start,end,lstart,lend,Ph,Scr);
  this->Next = ABlock;
  ABlock->Prev = this;
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
  Start = 0;
  End = 0;
  Evalue = 0.0;
  Length = 0;
  Strand = 0;
  NGaps = 0;
  Name = NULL;
  Match = NULL;
  Next = NULL;
}

// ---------------------------------------------------------------------
//  Construct from a char*...
// ---------------------------------------------------------------------
Hits :: Hits  (char* name, int length, char strand, int deb,int fin,int ldeb,int lfin,int Ph,int Scr,double Prob)
{
  Rejected = 0;
  Name = new char[strlen(name)+1];
  strcpy(Name,name);
  Strand = strand;
  Length = length;
  Start = deb;
  End = fin;
  Evalue = Prob;
  Match = new Block(deb,fin,ldeb,lfin,Ph,Scr);
  NGaps = 0;
  Next = NULL;
}

// ---------------------------------------------------------------------
//  Read a table of Hits from a file
// ---------------------------------------------------------------------
Hits* Hits::ReadFromFile(FILE* HitFile, int *NumHits) 
{
  char *HitId, *PHitId;
  int deb,fin,phase,HSPDeb,HSPFin,poids,read;
  double evalue,Pevalue;
  char A[128],B[128];
  Block *ThisBlock = NULL;
  Hits *OneHit = NULL, *ThisHit = NULL, *AllHit = NULL;

  Pevalue = -1.0;
  A[0] = B[0]= 0;
  HitId = A;
  PHitId = B;

  while ((read=fscanf(HitFile,"%d %d %d %lf %d %s %d %d\n",
		      &deb, &fin,&poids,&evalue,&phase,HitId,&HSPDeb,&HSPFin)) == 8)
    {
      if ((strcmp(HitId,PHitId) == 0) && 
	  (HSPDeb > ThisBlock->LEnd) && 
	  (deb > ThisBlock->End) &&
	  (evalue == Pevalue))
	// si HitId et PHitId sont égaux, alors il y a un Hit en cours
	// de meme nom on verifie que c'est bien compatible en terme
	// de position (sur l'est et le genomique) et en e-value.
	  {
	    ThisHit->NGaps++;
	    ThisHit->End = fin-1;
	    ThisBlock->AddBlockAfter(deb-1,fin-1,HSPDeb,HSPFin,phase,poids);
	    ThisBlock = ThisBlock->Next;
	  }
      else {
	(*NumHits)++;
	OneHit = new Hits(HitId,poids,(phase>0 ? 1 : -1),deb-1,fin-1,HSPDeb,HSPFin,phase,poids,evalue);
	ThisBlock = OneHit->Match;
	PHitId = HitId;
	Pevalue = evalue;

	if (AllHit == NULL) {
	  AllHit = OneHit;
	  ThisHit = OneHit;
	}
	else {
	  ThisHit->Next = OneHit;
	  ThisHit = OneHit;
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
  delete Name;
  delete Next;
}

