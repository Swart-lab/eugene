#ifndef  SENSOR_USER_H_INCLUDED
#define  SENSOR_USER_H_INCLUDED

#include "../../EuGene/Sensor.h"
#include "structure.h"
/*************************************************************
 **                       SensorUser                        **
 *************************************************************/
class SensorUser : public Sensor
{
 private:
  
  ptUTIL Signals;
  ptUTIL Contents;
  
 public:
  SensorUser  (int n, DNASeq *X);
  virtual ~SensorUser     ();
  virtual void Init       (DNASeq *);
  virtual void GiveInfo   (DNASeq *, int, DATA *);
  virtual void Plot       (DNASeq *);
  virtual void PostAnalyse(Prediction *);
};

extern "C" SensorUser * builder0( int n, DNASeq *X) {  return new SensorUser(n, X);}

#endif
