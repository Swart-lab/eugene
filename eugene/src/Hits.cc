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
Block ::  Block(int start, int end,int lstart,int lend)
{
  Start = start;
  End = end;
  LStart = lstart;
  LEnd = lend;
  Prev = Next = NULL;
}

// ---------------------------------------------------------------------
//  Insere apres
// ---------------------------------------------------------------------
void Block ::  AddBlockAfter(int start,int end,int lstart,int lend)
{
  Block * ABlock = new Block(start,end,lstart,lend);
  this->Next = ABlock;
  ABlock->Prev = this;
}
// ---------------------------------------------------------------------
//  Destroy this Block and all block after him
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
Hits :: Hits  (char* name, int length, char strand, int deb,int fin,int ldeb,int lfin)
{
  Rejected = 0;
  Name = new char[strlen(name)+1];
  strcpy(Name,name);
  Strand = strand;
  Length = length;
  Start = deb;
  Match = new Block(deb,fin,ldeb,lfin);
  NGaps = 0;
  Next = NULL;
}

// ---------------------------------------------------------------------
//  Read a table of Hits from a file
// ---------------------------------------------------------------------
Hits* Hits::ReadFromFile(FILE* ESTFile, int *NumEST) 
{
  char *EstId, *PEstId,*tmp;
  int deb,fin,brin,EstDeb,EstFin,poids;
  char A[128],B[128];
  Block *ThisBlock = NULL;
  Hits *OneEST = NULL,*ThisEST,*AllEST = NULL;

  A[0] = B[0]= 0;
  EstId = A;
  PEstId = B;

  while (fscanf(ESTFile,"%d %d %d %*s %d %s %d %d\n",
		&deb, &fin,&poids,&brin,EstId,&EstDeb,&EstFin) != EOF)
    {
      if ((strcmp(EstId,PEstId) == 0) && 
	  (EstDeb > ThisBlock->LEnd) && 
	  (deb > ThisBlock->End))
	// si EstId et PEstId sont égaux, alors il y a un EST en cours
	// de meme nom on verifie que c'est bien compatible en terme
	// de position (sur l'est et le genomique)
	  {
	    ThisEST->NGaps++;
	    ThisBlock->AddBlockAfter(deb-1,fin-1,EstDeb,EstFin);
	    ThisBlock = ThisBlock->Next;
	  }
      else {
	(*NumEST)++;
	OneEST = new Hits(EstId,poids,brin,deb-1,fin-1,EstDeb,EstFin);
	ThisBlock = OneEST->Match;
	tmp = PEstId;
	PEstId = EstId;
	EstId = tmp;

	if (AllEST == NULL) {
	  AllEST = OneEST;
	  ThisEST = OneEST;
	}
	else {
	  ThisEST->Next = OneEST;
	  ThisEST = OneEST;
	}
      }
    }
  return AllEST;
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

