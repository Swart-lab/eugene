#ifndef  SENSOR_EU_H_INCLUDED
#define  SENSOR_EU_H_INCLUDED

#include "../../EuGene/Sensor.h"

/*************************************************************
 **                     SensorEuStop                        **
 *************************************************************/
class SensorEuStop : public Sensor
{
 private:
  double   stopP;

 public:
  SensorEuStop  (int);
  virtual ~SensorEuStop   ();
  virtual void Init       (DNASeq *);
  virtual void GiveInfo   (DNASeq *X, int, DATA *);
  virtual void Plot       (DNASeq *X);
  virtual void PostAnalyse(Prediction *);
};

extern "C" SensorEuStop * builder0(int n) {  return new SensorEuStop(n); }

#endif
