#ifndef  SENSOR_EU_H_INCLUDED
#define  SENSOR_EU_H_INCLUDED

#include "Sensor.h"

/*************************************************************
 **                     SensorEuStop                        **
 *************************************************************/
class SensorEuStop : public Sensor
{
 private:
  REAL *Stop[2];

 public:
  SensorEuStop  ();
  ~SensorEuStop ();
  void Init     (DNASeq *);
  void GiveInfo (DNASeq *, int, DATA *);
};

#endif
