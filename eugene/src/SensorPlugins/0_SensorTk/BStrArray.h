// ---------------------------------------------------------------
//   T. Schiex
//
//     File:  BStrArray.h
//     Version:  1.0
//
//    Copyright (c) 2000 by Thomas Schiex All rights reserved.
//    Redistribution is not permitted without the express written
//    permission of the author.
//
//  Access to binary IMM model
// ---------------------------------------------------------------

#ifndef  BSTRARRAY_H_INCLUDED
#define  BSTRARRAY_H_INCLUDED

#include "../../EuGene/DNASeq.h"

class  BString_Array
{
 private:
  unsigned int  Max_Str_Len;
  int  Alphabet_Size;
  int  Num_Entries;
  unsigned short  * Val;
  
 public:
  int *Offset;
  BString_Array  ();
  BString_Array  (int, int);
  ~ BString_Array  ();
  int  Read  (FILE *);    
  int  String_To_Sub  (DNASeq  *, unsigned int, unsigned int);
  int  AntiString_To_Sub  (DNASeq  *, unsigned int, unsigned int);

  unsigned short operator [] (int i) {
    return  Val[i];
  }
  unsigned short &  operator ()  (DNASeq *, unsigned int, unsigned int);
};

#endif
