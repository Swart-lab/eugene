// ---------------------------------------------------------------
//   T. Schiex
//
//     File:  DNAseq.h
//     Version:  1.0
//
//    Copyright (c) 2000 by Thomas Schiex All rights reserved.
//    Redistribution is not permitted without the express written
//    permission of the authors.
//
//  Definitions for a class representing DNA seq with ambiguous data
//  and 6 possible unambiguous versions (6 phases possible completions)
// ---------------------------------------------------------------


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
  unsigned short Nuc2Code(char Ch);

 public:
  DNASeq ();
  DNASeq  (char* data);
  DNASeq (int size);
  DNASeq (FILE *);

  ~DNASeq ();

  int  SeqLen;
  char *Name;

  double Markov0[4];
  
  void  Print(FILE *);
  void  PrintTr(FILE*,int, int,signed char);

  unsigned char Degeneracy(int i, int sens);
  double IsStop(int i,int sens);
  double IsStart(int i,int sens);
  double IsEStart(int i,int sens);
  
  double Markov(int);
  double MarkovR(int);
  double GC_AT(int);
  
  void Transfer(int Pos, int Len, char *To, int mode);
  
  unsigned short Unambit (int pos);

  char operator [] (int i);
  char operator () (int i);
  unsigned short operator () (int i, int mode);
  char AA(int i, int mode);
  char nt(int i, int mode);
};
#endif
