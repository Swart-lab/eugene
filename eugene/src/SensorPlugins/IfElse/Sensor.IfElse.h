#ifndef  SENSOR_IFELSE_H_INCLUDED
#define  SENSOR_IFELSE_H_INCLUDED

#include "../Sensor.h"
#include "../Dll.h"

/*************************************************************
 **                     SensorIfElse                        **
 *************************************************************/
class SensorIfElse : public Sensor
{
 private:
  SensorLoader *sensorLIf;
  SensorLoader *sensorLElse;
  Sensor *sensorIf;
  Sensor *sensorElse;

 public:
  SensorIfElse  (int);
  virtual ~SensorIfElse   ();
  virtual void Init       (DNASeq *);
  virtual void GiveInfo   (DNASeq *, int, DATA *);
  virtual void Plot(DNASeq *);
  virtual void PostAnalyse(Prediction *);
};

extern "C" SensorIfElse * builder0(int n) {  return new SensorIfElse(n); }

#endif
