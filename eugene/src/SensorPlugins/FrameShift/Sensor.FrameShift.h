#ifndef  SENSOR_STARTCONST_INCLUDED
#define  SENSOR_STARTCONST_INCLUDED

#include "../Sensor.h"

/*************************************************************
 **                     SensorFrameShift                    **
 *************************************************************/
class SensorFrameShift : public Sensor
{
 private:
  // Probability of deletion/insertion
  REAL insProb;
  REAL delProb;

 public:
  SensorFrameShift  (int);
  virtual ~SensorFrameShift   ();
  virtual void Init       (DNASeq *);
  virtual void GiveInfo   (DNASeq *, int, DATA *);
  virtual void Plot       (DNASeq *);
  virtual void PostAnalyse(Prediction *);
};

extern "C" SensorFrameShift * builder0(int n) {  return new SensorFrameShift(n); }

#endif
