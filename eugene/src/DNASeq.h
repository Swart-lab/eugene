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
// File:     DNASeq.cc
// Contents: Definitions for a class representing DNA seq with ambiguous data
//  and 6 possible unambiguous versions (6 phases possible completions)
// ------------------------------------------------------------------


#ifndef  DNASEQ_H_INCLUDED
#define  DNASEQ_H_INCLUDED

#include <stdio.h>

const unsigned short MASKSEQ = 0x000F;

const unsigned short BitT = 3;
const unsigned short BitG = 2;
const unsigned short BitC = 1;
const unsigned short BitA = 0;


const unsigned short CodeT = (1 << BitT);
const unsigned short CodeG = (1 << BitG);
const unsigned short CodeC = (1 << BitC);
const unsigned short CodeA = (1 << BitA);

class  DNASeq
{
 private:
  int  Size;
  unsigned short *Sequence;

  void UpdateMarkov();

 public:
  DNASeq ();
  DNASeq (int size);
  DNASeq (char* filename);

  ~DNASeq ();

  int  SeqLen;
  char *Name;

  double Markov0[4];
  
  void  Print(FILE *);
  void  PrintTr(FILE*,int, int,signed char);

  unsigned char Degeneracy(int i, int sens);

  const int static StopAfterT  = 1;
  const int static StopAfterTG = 2;
  const int static StopAfterTA = 4;
  const int static StopAfterG = 8;
  const int static StopAfterA = 16;
  const int static StopAfterAG = 32;

  int    IsStartStop(int i);
  int    IsStopStop(int i);
  double IsStop(int i,int sens);
  double IsStart(int i,int sens);
  double IsEStart(int i,int sens);
  
  double Markov(int);
  double MarkovR(int);
  double GC_AT(int);

  unsigned short Nuc2Code(char Ch);

  void Transfer(int Pos, int Len, char *To, int mode);
  
  unsigned short Unambit (int pos);

  char operator [] (int i);
  char operator () (int i);
  unsigned short operator () (int i, int mode);
  char AA(int i, int mode);
  char nt(int i, int mode);
};
#endif
