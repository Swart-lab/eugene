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
// File:     Prediction.cc
// Contents: class Prediction
// ------------------------------------------------------------------

#include<iostream>

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
//  add.
// ------------------------
void Prediction :: poptill (int pos)
{

  while (vPos.back() < pos) {
    vPos.pop_back();
    vState.pop_back();
    nb--;
  }
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
	     1+(getNextState(i)==InterGen)*3);
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
  int stUtr = 0;
  int utr53 = 0;

  // ATTENTION aux utrs en debut de prediction.
  while(vState[i] >= InterGen) {
    if (vState[i] == UTR5F) {
      if (stUtr == -1  ||  utr53 == 3) break;
      stUtr = 1;
      utr53 = 5;
    }
    if (vState[i] == UTR5R) {
      if (stUtr == 1  ||  utr53 == 3) break;
      stUtr = -1;
      utr53 = 5;
    }
    if (vState[i] == UTR3F) {
      if (stUtr == -1  ||  utr53 == 5) break;
      stUtr = 1;
      utr53 = 3;
    }
    if (vState[i] == UTR3R) {
      if (stUtr == 1  ||  utr53 == 5) break;
      stUtr = -1;
      utr53 = 3;
    }
    i--;
  }
  for(int j=1; j<geneNumber; j++) {
    while(vState[i] <= InterGen)
      i--;
    while(vState[i] >= InterGen)
      i--;
  }
  while(vState[i] <= InterGen && i != -1) {
    if(vState[i] <= TermR3)
      nb++;
    i--;
  }
  return nb;
}

// ------------------------
//  lenCDS.
// ------------------------
int Prediction :: lenCDS (int geneNumber)
{
  int i     = (int)vPos.size()-1;
  int len   = 0;
  int stUtr = 0;
  int utr53 = 0;
  
  // ATTENTION aux utrs en debut de prediction.
  while(vState[i] >= InterGen) {
    if (vState[i] == UTR5F) {
      if (stUtr == -1  ||  utr53 == 3) break;
      stUtr = 1;
      utr53 = 5;
    }
    if (vState[i] == UTR5R) {
      if (stUtr == 1  ||  utr53 == 3) break;
      stUtr = -1;
      utr53 = 5;
    }
    if (vState[i] == UTR3F) {
      if (stUtr == -1  ||  utr53 == 5) break;
      stUtr = 1;
      utr53 = 3;
    }
    if (vState[i] == UTR3R) {
      if (stUtr == 1  ||  utr53 == 5) break;
      stUtr = -1;
      utr53 = 3;
    }
    i--;
  }
  for(int j=1; j<geneNumber; j++) {
    while(vState[i] <= InterGen)
      i--;
    while(vState[i] >= InterGen)
      i--;
  }
  while(vState[i] <= InterGen && i != -1) {
    if(vState[i] <= TermR3)
      len += (vPos[i] - vPos[i+1]-1) + 1;
    i--;
  }
  return len;
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
  int state  = getStateForPos (p);
  int nState = getStateForPos (p+1);
  
  if(1 <= State2Phase[nState] && State2Phase[nState] <= 3  &&  state >= InterGen)
    return "True";
  if(-3 <= State2Phase[state] && State2Phase[state] <= -1  && 
     (nState >= InterGen || p+1 >= vPos[0]))
    return "True";
  return "False";
}

// ------------------------
//  isStop.
// ------------------------
char* Prediction :: isStop(int p)
{
  int state  = getStateForPos (p);
  int nState = getStateForPos (p+1);
  
  if(state != -1  &&  1 <= State2Phase[state] && State2Phase[state] <= 3  &&  
     (nState >= InterGen || p+1 >= vPos[0]))
    return "True";
  if(-3 <= State2Phase[nState] && State2Phase[nState] <= -1 &&  state >= InterGen)
    return "True";
  return "False";
}

// ------------------------
//  isDon.
// ------------------------
char* Prediction :: isDon(int p)
{
  int state  = getStateForPos (p);
  int nState = getStateForPos (p+1);
  
  if(1 <= State2Phase[state] && State2Phase[state] <= 3 &&  nState == IntronF1)
    return "True";
  if(-3 <= State2Phase[nState] && State2Phase[nState] <= -1 &&  state == IntronR1)
    return "True";
  return "False";
}

// ------------------------
//  isAcc.
// ------------------------
char* Prediction :: isAcc(int p)
{
  int pState = getStateForPos (p);
  int state  = getStateForPos (p+1);
  
  if(1 <= State2Phase[state] && State2Phase[state] <= 3 &&  pState == IntronF1)
    return "True";
  if(-3 <= State2Phase[pState] && State2Phase[pState] <= -1  &&  state == IntronR1)
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
//         is InterGen and the state after the last exon 
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
  int state, nState, pState;

  if (sig_type == DATA::Start) {
    state  = getStateForPos (pos);
    nState = getStateForPos (pos+1);
    if (strand == '+') {
      if ((1 <= State2Phase[nState] && State2Phase[nState] <= 3) && state == InterGen) 
	is_state = true;
    } else if (strand == '-') {
      if ((-3 <= State2Phase[state] && State2Phase[state] <= -1)  &&  
	  (nState == UTR5R ||  pos == vPos[0])) 
	is_state = true;
    } else bad_strand = true;

  } else if (sig_type == DATA::Stop) {
    state  = getStateForPos (pos);
    nState = getStateForPos (pos+1);
    if (strand == '+') {
      if( (1 <= State2Phase[state] && State2Phase[state] <= 3) &&  
	  (nState == UTR3F || pos == vPos[0])) is_state = true;
    } else if (strand == '-') {
      if( (-3 <= State2Phase[nState] && State2Phase[nState] <= -1)  &&  
	  (state == UTR3R || state == InterGen) ) is_state = true;
    } else bad_strand = true;

  } else if (sig_type == DATA::Acc) {
    pState  = getStateForPos (pos);
    state = getStateForPos (pos+1);
    if (strand == '+') {
      if ((1 <= State2Phase[state] && State2Phase[state]<= 3) &&  
	  (4 <= State2Frame[pState] && State2Frame[pState] <= 6)) is_state = true;
    } else if (strand == '-') {
      if( (-3 <= State2Phase[pState] && State2Phase[pState] <= -1)  &&  
	  (-6 <= State2Frame[state] && State2Frame[state] <= -4)) is_state = true;
    } else bad_strand = true;

  } else if (sig_type == DATA::Don) {
    state  = getStateForPos (pos);
    nState = getStateForPos (pos+1);
    if (strand == '+') {
      if ((1 <= State2Phase[state] && State2Phase[state] <= 3) &&  
	  (4 <= State2Frame[nState] && State2Frame[nState] <= 6)) is_state = true;
    } else if (strand == '-') {
      if( (-3 <= State2Phase[nState] && State2Phase[nState] <= -1) &&  
	  (-6 <= State2Frame[state] && State2Frame[state] <= -4)) is_state = true;
    } else bad_strand = true;

  } else 
    {std::cerr<<"ERROR: bad state "<<sig_type<<" given in argument in Prediction::IsState.\n"; exit(2);}

  if (bad_strand)
    {std::cerr<<"ERROR: bad strand "<<strand<<"  given in argument in Prediction::IsState.\n"; exit(2);}

  return is_state;
}


