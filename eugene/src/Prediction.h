#ifndef  PREDICTION_H_INCLUDED
#define  PREDICTION_H_INCLUDED

#include <vector>
#include "Const.h"
extern "C"{
#include "../GDIF/gdIF.h"
}

/*************************************************************
 **                        Prediction                       **
 *************************************************************/
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
  void add          (int, signed char);
  void setPos       (int, int);
  char getNextState (int);
  char getState     (int i) { return vState[i]; }
  int  getPos       (int i) { return vPos[i];   }
  int  size         ()      { return nb;        } 
  void plotPred     ();
  void resetPred    ();
  int  nbExon       (int);
};

#endif
