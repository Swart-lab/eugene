#ifndef  SENSOR_EU_H_INCLUDED
#define  SENSOR_EU_H_INCLUDED

#include "../Sensor.h"

/*************************************************************
 **                     SensorEuStop                        **
 *************************************************************/
class SensorEuStop : public Sensor
{
 private:
  REAL *Stop[2];

 public:
  SensorEuStop  (int);
  virtual ~SensorEuStop ();
  virtual void Init     (DNASeq *);
  virtual void GiveInfo (DNASeq *, int, DATA *);
};

extern "C" SensorEuStop * builder0(int n) {  return new SensorEuStop(n); }

#endif
