#include "Prediction.h"

// ------------------------
//  Default constructor.
// ------------------------
Prediction :: Prediction ()
{
  index = 0;
  nb = 0;
}

// ------------------------
//  Default destructor.
// ------------------------
Prediction :: ~Prediction ()
{
  vPos.clear();
  vState.clear();
}

// ------------------------
//  add.
// ------------------------
void Prediction :: add (int pos, signed char state)
{
  vPos.push_back  ( pos   );
  vState.push_back( state );
  nb++;
}

// ------------------------
//  print prediction.
// ------------------------
void Prediction :: print ()
{
  for(int i=(int)vPos.size()-1; i!=-1; i--)
    printf("State:%2d\tPos:%d\n",vState[i],vPos[i]);
}

// ------------------------
//  getNextState.
// ------------------------
char Prediction :: getNextState (int pos)
{
  if(vPos[index + nb-1] < pos)
    index--;
  return vState[index + nb-1];
}

// ------------------------
//  getStateForPos.
// ------------------------
char Prediction :: getStateForPos(int pos)
{
  int i = vPos.size()-1;
  while(i!=-1 && vPos[i] < pos)
    i--;
  if(i!=-1) return vState[i];
  else      return -1;
}

// ------------------------
//  setPos.
// ------------------------
void Prediction :: setPos (int i, int newPos)
{ vPos[i] = newPos; }

// ------------------------
//  plotPred.
// ------------------------
void Prediction :: plotPred ()
{
  const int PredWidth = 2;
  const short int State2Phase[18] = {1,2,3,-1,-2,-3,4,4,4,-4,-4,-4,0,0,0,0,0,0};
  int SeqLen = vPos[0] - 1;
  index = 0;
  
  for (int i=0; i<SeqLen; i++)
    PlotBarI(i, State2Phase[getNextState(i)], 0.4, PredWidth,
	     1+((getNextState(i)==12) || (getNextState(i)==17))*3);
}

// ------------------------
//  resetPred.
// ------------------------
void Prediction :: resetPred ()
{
  index = 0;
  nb = 0;
  vPos.clear();
  vState.clear();
}

// ------------------------
//  nbExon.
// ------------------------
int Prediction :: nbExon (int geneNumber)
{
  int i = (int)vPos.size()-1;
  int nb = 0;
  
  while(vState[i] >= InterGen5)
      i--;
  
  for(int j=1; j<geneNumber; j++) {
    while(vState[i] <= InterGen5)
      i--;
    while(vState[i] >= InterGen5)
      i--;
  }
  while(vState[i] <= InterGen5 && i != -1) {
    if(vState[i] <= ExonR3)
      nb++;
    i--;
  }
  return nb;
}
