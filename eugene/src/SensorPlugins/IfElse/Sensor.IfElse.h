#ifndef  SENSOR_IFELSE_H_INCLUDED
#define  SENSOR_IFELSE_H_INCLUDED

#include "../../EuGene/Sensor.h"

/*************************************************************
 **                     SensorIfElse                        **
 *************************************************************/
class SensorIfElse : public Sensor
{
 private:
  Sensor *sensorIf;
  Sensor *sensorElse;

 public:
  SensorIfElse  (int n, DNASeq *X);
  virtual void Init       (DNASeq *);
  virtual void GiveInfo   (DNASeq *, int, DATA *);
  virtual void Plot(DNASeq *);
  virtual void PostAnalyse(Prediction *);
};

extern "C" SensorIfElse * builder0(int n, DNASeq *X) {  return new SensorIfElse(n, X); }

#endif
