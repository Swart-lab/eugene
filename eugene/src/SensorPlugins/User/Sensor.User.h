#ifndef  SENSOR_USER_H_INCLUDED
#define  SENSOR_USER_H_INCLUDED

#include "../Sensor.h"

/*************************************************************
 **                       SensorUser                        **
 *************************************************************/
class SensorUser : public Sensor
{
 private:
  
 public:
  SensorUser  (int);
  virtual ~SensorUser     ();
  virtual void Init       (DNASeq *);
  virtual void ResetIter  ();
  virtual void GiveInfo   (DNASeq *, int, DATA *);
  virtual void GiveInfoAt (DNASeq *, int, DATA *);
  virtual void Plot       (DNASeq *);
  virtual void PostAnalyse(Prediction *);
};

extern "C" SensorUser * builder0( int n ) {  return new SensorUser(n);}

#endif
