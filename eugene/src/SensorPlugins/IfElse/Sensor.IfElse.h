#ifndef  SENSOR_IFELSE_H_INCLUDED
#define  SENSOR_IFELSE_H_INCLUDED

#include "../../EuGene/Sensor.h"
#include "../../EuGene/Dll.h"

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
  SensorIfElse  (int n, DNASeq *X);
  virtual ~SensorIfElse   ();
  virtual void Init       (DNASeq *);
  virtual void GiveInfo   (DNASeq *, int, DATA *);
  virtual void Plot(DNASeq *);
  virtual void PostAnalyse(Prediction *);
};

extern "C" SensorIfElse * builder0(int n, DNASeq *X) {  return new SensorIfElse(n, X); }

#endif
