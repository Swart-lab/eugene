#include "Prediction.h"

#include<iostream>

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

 
// ------------------------
// IsState: Is a nucleotid in a given state ? 
// BE CAREFUL: 
// this method could be used only for a prediction with one complete gene.
//         the prediction could be the representation of an external gff annotation 
//         that could not specify the UTR. In this case, the state from 0 to the first exon
//         in the annotation (the first element of vPos and vState)
//         is InterGen5 and the state after the last exon 
//         in the annotation (the last element of vPos and vState)
//         is not set (getStateForPos returns -1).
//         To identify Start Reverse and Stop Forward, an other condition is used
//         pos == vPos[0] that compares the position with 
//                             in Forward the end of the last exon of the gene (vPos[0])
//                             in Reverse the start of the first exon of the gene (vPos[0])
// ------------------------
bool Prediction :: IsState (DATA::SigType sig_type, int pos, char strand)
{
  bool is_state = false;
  bool bad_strand = false;
  char state, nState, pState;

  if (sig_type == DATA::Start) {
    state  = getStateForPos (pos);
    nState = getStateForPos (pos+1);
    if (strand == '+') {
      if( (nState == ExonF1 || nState == ExonF2 || nState == ExonF3 ) &&  
	  state == InterGen5) is_state = true;
    } else if (strand == '-') {
      if( (state == ExonR1 || state == ExonR2 || state == ExonR3)  &&  
	  (nState == UTR5R ||  pos == vPos[0]) ) is_state = true;
    } else bad_strand = true;

  } else if (sig_type == DATA::Stop) {
    state  = getStateForPos (pos);
    nState = getStateForPos (pos+1);
    if (strand == '+') {
      if( (state == ExonF1 || state == ExonF2  || state == ExonF3) &&  
	  (nState == UTR3F || pos == vPos[0])) is_state = true;
    } else if (strand == '-') {
      if( (nState == ExonR1 || nState == ExonR2 || nState == ExonR3)  &&  
	  (state == UTR3R || state == InterGen5) ) is_state = true;
    } else bad_strand = true;

  } else if (sig_type == DATA::Acc) {
    pState  = getStateForPos (pos);
    state = getStateForPos (pos+1);
    if (strand == '+') {
      if( (state == ExonF1 || state == ExonF2 || state == ExonF3)  &&  
	  (pState == IntronF1 || pState == IntronF2 || pState == IntronF3) ) is_state = true;
    } else if (strand == '-') {
      if( (pState == ExonR1 || pState == ExonR2 || pState == ExonR3)  &&  
	  (state == IntronR1 || state == IntronR2 || state == IntronR3) ) is_state = true;
    } else bad_strand = true;

  } else if (sig_type == DATA::Don) {
    state  = getStateForPos (pos);
    nState = getStateForPos (pos+1);
    if (strand == '+') {
      if( (state == ExonF1 || state == ExonF2 || state == ExonF3) &&  
	  (nState == IntronF1 || nState == IntronF2 || nState == IntronF3) ) is_state = true;
    } else if (strand == '-') {
      if( (nState == ExonR1 || nState == ExonR2 || nState == ExonR3) &&  
	  (state == IntronR1 || state == IntronR2 || state == IntronR3) ) is_state = true;
    } else bad_strand = true;

  } else 
    {std::cerr<<"ERROR: bad state "<<sig_type<<" given in argument in Prediction::IsState.\n"; exit(2);}

  if (bad_strand)
    {std::cerr<<"ERROR: bad strand "<<strand<<"  given in argument in Prediction::IsState.\n"; exit(2);}

  return is_state;
}


