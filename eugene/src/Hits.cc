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
  Length = 0;
  Strand = 0;
  Name = NULL;
  Match = NULL;
  Next = NULL;
}

// ---------------------------------------------------------------------
//  Construct from a char*...
// ---------------------------------------------------------------------
Hits :: Hits  (char* name, int length, char strand, int deb,int fin,int ldeb,int lfin)
{

  Name = new char[strlen(name)+1];
  strcpy(Name,name);
  Strand = strand;
  Length = length;
  Match = new Block(deb,fin,ldeb,lfin);
  Next = NULL;
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

