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

extern "C"{
#include "../GDIF/gdIF.h"
}

/*************************************************************
 **                        Sensor                           **
 *************************************************************/
class Sensor
{
 private:
  
 public:
  TYPE_SENSOR type;

  Sensor  ();
  virtual ~Sensor ();
  virtual void Init     (DNASeq *) = 0;
  virtual void GiveInfo (DNASeq *, int, DATA *) = 0;
  void CheckStart   (DNASeq *, REAL **);
  void CheckSplices (DNASeq *, REAL **, REAL **);
};

class SensorFactory
{
 public:
  SensorFactory() { }
  virtual ~SensorFactory() { }
  virtual Sensor * CreateSensor() = 0;
};
#endif