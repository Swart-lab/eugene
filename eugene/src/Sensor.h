#ifndef  SENSOR_H_INCLUDED
#define  SENSOR_H_INCLUDED

#include <string.h>
#ifdef HAVE_STRINGS_H
#include <strings.h>
#endif

#include "Param.h"
#include "DNASeq.h"
#include "Plot.h"
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
  virtual void Init     (DNASeq *) = 0;
  virtual void GiveInfo (DNASeq *, int, DATA *) = 0;
  void CheckStart   (DNASeq *, REAL **);
  void CheckSplices (DNASeq *, REAL **, REAL **);
  int GetNumber() { return instanceNumber; }
};

#endif

