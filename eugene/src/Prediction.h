#ifndef  PREDICTION_H_INCLUDED
#define  PREDICTION_H_INCLUDED

#include <cstdio>
#include <vector>
#include <algorithm>
#include "Const.h"
extern "C"{
#include "../GDIF/gdIF.h"
}

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
};

// Les etats
enum Tracks { 
  ExonF1 = 0, ExonF2 = 1, ExonF3 = 2,
  ExonR1 = 3, ExonR2 = 4, ExonR3 = 5,
  IntronF1 = 6,IntronF2 = 7, IntronF3 = 8,
  IntronR1 = 9, IntronR2 = 10, IntronR3 = 11,
  InterGen5 = 12, InterGen3 = 17,
  UTR5F = 13, UTR5R = 15,
  UTR3F = 14, UTR3R = 16};
const short int State2Phase[18] = {1,2,3,-1,-2,-3,4,4,4,-4,-4,-4,0,0,0,0,0,0};

#endif
