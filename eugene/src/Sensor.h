#ifndef  SENSOR_H_INCLUDED
#define  SENSOR_H_INCLUDED

#include <vector>
#include <algorithm>
#include <string.h>
#ifdef HAVE_STRINGS_H
#include <strings.h>
#endif

#include "Param.h"
#include "DNASeq.h"
#include "Prediction.h"
#include "System.h"
#include "Const.h"

extern "C" {
#include "../GDIF/gdIF.h"
}

/*************************************************************
 **                        Sensor                           **
 *************************************************************/
class Sensor
{
 private:
  int instanceNumber;
  
 public:
  TYPE_SENSOR type;
  
  Sensor  (int);
  virtual ~Sensor ();
  virtual void Init       (DNASeq *) = 0;
  virtual void ResetIter  () = 0;
  virtual void GiveInfo   (DNASeq *, int, DATA *) = 0;
  virtual void GiveInfoAt (DNASeq *, int, DATA *) = 0;
  virtual void Plot       (DNASeq *) = 0;
  virtual void PostAnalyse(Prediction *) = 0;
  void CheckStart   (DNASeq *, vector<int>, vector<int>);
  void CheckSplices (DNASeq *, vector<int>, vector<int>, vector<int>, vector<int>);
  int  GetNumber    () { return instanceNumber; }
};

#endif

