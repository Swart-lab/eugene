#include "Sensor.h"

/*************************************************************
 **                        Sensor                           **
 *************************************************************/
// ----------------------
//  Default constructor.
// ----------------------
Sensor :: Sensor (int n)
{
  instanceNumber = n;
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
void Sensor :: CheckStart(DNASeq *X, std::vector<int> vPosF, std::vector<int> vPosR)
{
  for (int i = 0; i<(int)vPosF.size(); i++)
    if (!X->IsEStart(vPosF[i],1))
      fprintf(stderr,"WARNING: Non ATG start predicted at %d on + strand!\n", vPosF[i]+1);
  
  for (int i = 0; i<(int)vPosR.size(); i++)
    if (!X->IsEStart(vPosR[i]-1,-1))
      fprintf(stderr,"WARNING: Non ATG start predicted at %d on - strand!\n", vPosR[i]+1);
}

// -----------------------------------------------
//  Sanity check; AG/GT or GC  splice site check.
// -----------------------------------------------
void Sensor :: CheckSplices (DNASeq *X,
			     std::vector<int> vPosAccF, std::vector<int> vPosDonF,
			     std::vector<int> vPosAccR, std::vector<int> vPosDonR)
{
  for (int i = 0; i<(int)vPosAccF.size(); i++)
    if(((*X)[vPosAccF[i]-2] != 'a') || ((*X)[vPosAccF[i]-1] != 'g'))
      fprintf(stderr,"WARNING: Non AG (%c%c) acceptor at %d (+ strand) !\n", 
	      (*X)[vPosAccF[i]-2], (*X)[vPosAccF[i]-1], vPosAccF[i]);
  
  for (int i = 0; i<(int)vPosAccR.size(); i++)
    if(((*X)(vPosAccR[i]) != 'g') || ((*X)(vPosAccR[i]+1) != 'a'))
      fprintf(stderr,"WARNING: Non AG (%c%c) acceptor at %d (- strand) !\n",
	      (*X)(vPosAccR[i]), (*X)(vPosAccR[i]+1),vPosAccR[i]);
  
  for (int i = 0; i<(int)vPosDonF.size(); i++)
    if(((*X)[vPosDonF[i]] != 'g') || 
       (((*X)[vPosDonF[i]+1] != 't') && ((*X)[vPosDonF[i]+1] != 'c')))
      fprintf(stderr,"WARNING: Non GT/GC (%c%c) donor at %d (+ strand) !\n", 
	      (*X)[vPosDonF[i]],(*X)[vPosDonF[i]+1], vPosDonF[i]);
  
  for (int i = 0; i<(int)vPosDonR.size(); i++)
    if((((*X)(vPosDonR[i]-2) != 't') &&
	((*X)(vPosDonR[i]-2) != 'c')) || ((*X)(vPosDonR[i]-1) != 'g'))
      fprintf(stderr,"WARNING: Non GT/GC (%c%c) donor at %d (- strand) !\n", 
	      (*X)(vPosDonR[i]-2),(*X)(vPosDonR[i]-2),vPosDonR[i]);
}
