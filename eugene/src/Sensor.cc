#include "Sensor.h"

/*************************************************************
 **                        Sensor                           **
 *************************************************************/
// ----------------------
//  Default constructor.
// ----------------------
Sensor :: Sensor ()
{
  type = Type_Unknown;
}

// ----------------------
//  Default destructor.
// ----------------------
Sensor :: ~Sensor ()
{
}

// ---------------------------------------
//  Sanity check : check that ATG occurs.
// ---------------------------------------
void Sensor :: CheckStart(DNASeq *X, REAL **Start)
{
  for (int pos = 0; pos <= X->SeqLen; pos++) {
    if ((Start[0][pos] != 0.0) && (!X->IsEStart(pos,1)))
      fprintf(stderr,"WARNING: Non ATG start predicted at %d on + strand!\n",
	      pos+1);
    
    if ((Start[1][pos] != 0.0) && (!X->IsEStart(pos-1,-1)))
      fprintf(stderr,"WARNING: Non ATG start predicted at %d on - strand!\n",
	      pos+1);
  }
}

// -----------------------------------------------
//  Sanity check; AG/GT or GC  splice site check.
// -----------------------------------------------
void Sensor :: CheckSplices (DNASeq *X, REAL **Acc, REAL **Don)
{
  for (int pos = 0; pos <= X->SeqLen; pos++) {
    if (Acc[0][pos] != 0.0) 
      if(((*X)[pos-2] != 'a') || ((*X)[pos-1] != 'g'))
	fprintf(stderr,"WARNING: Non AG acceptor at %d (+ strand) !\n",pos);
    
    if (Acc[1][pos] != 0.0)
      if(((*X)(pos) != 'g') || ((*X)(pos+1) != 'a'))
	fprintf(stderr,"WARNING: Non AG acceptor at %d (- strand) !\n",pos);
    
    if (Don[0][pos] != 0.0)
      if(((*X)[pos] != 'g') || (((*X)[pos+1] != 't') && ((*X)[pos+1] != 'c')))
	fprintf(stderr,"WARNING: Non GT/GC donor at %d (+ strand) !\n",pos);
    
    if (Don[1][pos] != 0.0)
      if((((*X)(pos-2) != 't') && ((*X)(pos-2) != 'c')) || ((*X)(pos-1)!='g'))
	fprintf(stderr,"WARNING: Non GT/GC donor at %d (- strand) !\n",pos);
  }
}
