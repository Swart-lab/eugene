#ifndef  PREDICTION_H_INCLUDED
#define  PREDICTION_H_INCLUDED

#include <cstdio>
#include <vector>
#include <algorithm>

#include "Const.h"
extern "C"{
#include "../GDIF/gdIF.h"
}

#include "SensorIF.h"

/*************************************************************
 **                        Prediction                       **
 ************************************************************/
class Prediction
{
 private:
  int index;
  int nb;
  std::vector <int> vPos;
  std::vector <signed char> vState;
  
 public:
  Prediction  ();
  ~Prediction ();

  double OptimalPath;

  void  add           (int, signed char);
  void  print         ();
  void  setPos        (int, int);
  char  getNextState  (int); 
  char  getStateForPos(int);
  char  getState      (int i) { return vState[i]; }
  int   getPos        (int i) { return vPos[i];   }
  int   size          ()      { return nb;        } 
  void  plotPred      ();
  void  resetPred     ();
  int   nbExon        (int);
  void  reversePred   ();
  char* isStart       (int);
  char* isStop        (int);
  char* isDon         (int);
  char* isAcc         (int);
  bool IsState (DATA::SigType sig_type, int pos, char strand);
};

// Les etats
enum Tracks {
  InitF1 = 0, InitF2 = 1, InitF3 = 2,
  InitR1 = 3, InitR2 = 4, InitR3 = 5,
  SnglF1 = 6, SnglF2 = 7, SnglF3 = 8,
  SnglR1 = 9, SnglR2 = 10, SnglR3 = 11,
  IntrF1 = 12, IntrF2 = 13, IntrF3 = 14,
  IntrR1 = 15, IntrR2 = 16, IntrR3 = 17,
  TermF1 = 18, TermF2 = 19, TermF3 = 20,
  TermR1 = 21, TermR2 = 22, TermR3 = 23,
  IntronF1 = 24,IntronF2 = 25, IntronF3 = 26,
  IntronR1 = 27, IntronR2 = 28, IntronR3 = 29,
  InterGen = 30, 
  UTR5F = 31, UTR3F = 32, 
  UTR5R = 33,UTR3R = 34, 
  IntronU5F = 35, IntronU5R = 36, 
  IntronU3F = 37, IntronU3R = 38, 
  NbTracks = 39};

const short int State2Phase[NbTracks] = {1,2,3,-1,-2,-3,
					 1,2,3,-1,-2,-3,
					 1,2,3,-1,-2,-3,
					 1,2,3,-1,-2,-3,
					 4,4,4,
					 -4,-4,-4,
					 0,0,0,0,0,
					 4,-4,4,-4};

const short int State2Frame[NbTracks] = {1,2,3,-1,-2,-3,
					 1,2,3,-1,-2,-3,
					 1,2,3,-1,-2,-3,
					 1,2,3,-1,-2,-3,
					 4,5,6,
					 -4,-5,-6,
					 0,0,0,0,0,
					 4,-4,4,-4};

inline int PhaseAdapt(char p) {return State2Frame[(int)p];}

#endif
