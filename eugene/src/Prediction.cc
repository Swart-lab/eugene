#include "Prediction.h"

// ------------------------
//  Default constructor.
// ------------------------
Prediction :: Prediction ()
{
  index = 0;
  nb = 0;
  OptimalPath = 0;
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
  int SeqLen = vPos[0] - 1;
  index = 0;

  for (int i=0; i<SeqLen; i++)
    PlotBarI(i, State2Phase[(int)getNextState(i)], 0.4, PredWidth,
	     1+((getNextState(i)==12) || (getNextState(i)==17))*3);
}

// ------------------------
//  resetPred.
// ------------------------
void Prediction :: resetPred ()
{
  index = 0;
  nb = 0;
  OptimalPath = 0;
  vPos.clear();
  vState.clear();
}

// ------------------------
//  nbExon.
// ------------------------
int Prediction :: nbExon (int geneNumber)
{
  int i     = (int)vPos.size()-1;
  int nb    = 0;
  int nbUtr = 0;    // Pb UTR :
                    //  - seq.1.1.0 Utr5 - ...  Gene 1
                    //  - seq.1.2.0 Utr5 + ...  Gene 2
  while(vState[i] >= InterGen5) {
    if (vState[i] >= UTR5F  &&  vState[i] <= UTR3R)
      nbUtr++;
    if (nbUtr > 1)
      break;
    i--;
  }

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

// ----------------
//  reverse.
// ----------------
void Prediction :: reversePred()
{
  reverse(vState.begin(), vState.end());
  reverse(vPos.begin(),   vPos.end());
}

// ------------------------
//  isStart.
// ------------------------
char* Prediction :: isStart(int p)
{
  char state  = getStateForPos (p);
  char nState = getStateForPos (p+1);
  
  if(nState <= ExonF3  &&  state >= InterGen5)
    return "True";
  if(state >= ExonR1  &&  state <= ExonR3  && 
     (nState >= InterGen5 || p+1 >= vPos[0]))
    return "True";
  return "False";
}

// ------------------------
//  isStop.
// ------------------------
char* Prediction :: isStop(int p)
{
  char state  = getStateForPos (p);
  char nState = getStateForPos (p+1);
  
  if(state != -1  &&  state <= ExonF3  &&  
     (nState >= InterGen5 || p+1 >= vPos[0]))
    return "True";
  if(nState >= ExonR1  &&  nState <= ExonR3  &&  state >= InterGen5)
    return "True";
  return "False";
}

// ------------------------
//  isDon.
// ------------------------
char* Prediction :: isDon(int p)
{
  char state  = getStateForPos (p);
  char nState = getStateForPos (p+1);
  
  if(state <= ExonF3  &&  nState == IntronF1)
    return "True";
  if(nState >= ExonR1  &&  nState <= ExonR3  &&  state == IntronR1)
    return "True";
  return "False";
}

// ------------------------
//  isAcc.
// ------------------------
char* Prediction :: isAcc(int p)
{
  char pState = getStateForPos (p);
  char state  = getStateForPos (p+1);
  
  if(state <= ExonF3  &&  pState == IntronF1)
    return "True";
  if(pState >= ExonR1  &&  pState <= ExonR3  &&  state == IntronR1)
    return "True";
  return "False";
}
