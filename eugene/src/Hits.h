/*****************************************************************************/
/*             Copyright (c) 2004 by INRA. All rights reserved.              */
/*                 Redistribution is not permitted  without                  */
/*                 the express written permission of  INRA.                  */
/*                   Mail : eugene@ossau.toulouse.inra.fr                    */
/*---------------------------------------------------------------------------*/
/* File         : EuGeneTk/EuGene/Hits.h                                     */
/* Description  : Definitions for a class representing alignements with      */
/*                genomic seq.                                               */
/* Authors      : P.Bardou, S.Foissac, M.J.Cros, A.Moisan, T.Schiex          */
/* History      : July 2004                                                  */
/*****************************************************************************/

#ifndef  HITS_H_INCLUDED
#define  HITS_H_INCLUDED
#include <stdio.h>

class Block
{
 private:
  
 public:
  Block ();
  Block(int Start, int End, int LStart, int LEnd, int Ph, int Scr);
  ~Block ();
  
  void AddBlockAfter(int Start,int End,int LStart,int LEnd,int Ph,int Scr);
  
  int Start;
  int End;
  int Phase;
  int Score;
  int LStart;
  int LEnd;
  Block *Prev,*Next;
};

class  Hits
{
 private:

 public:
  Hits ();
  Hits (char* name, int length, char strand, int deb, int fin, int ldeb,
	int lfin, int Ph, int Scr, double prob, int level, int sup);
  Hits* Hits::ReadFromFile(FILE* HitFile, int *NumSeqs, int level, int margin);
  ~Hits ();

  char   *Name;
  char   Strand;
  int    Length;
  int    NGaps;
  int    Start;
  int    End;
  double Evalue;
  char   Rejected;
  int    Level;         // For blastx level 0-9 
  int    Support;

  Block *Match;
  Hits  *Next;
};
#endif
