#ifndef  SENSOR_STARTCONST_INCLUDED
#define  SENSOR_STARTCONST_INCLUDED

#include "../../EuGene/Sensor.h"

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
  SensorFrameShift  (int n, DNASeq *X);
  virtual ~SensorFrameShift   ();
  virtual void Init       (DNASeq *X);
  virtual void GiveInfo   (DNASeq *X, int, DATA *);
  virtual void Plot       (DNASeq *X);
  virtual void PostAnalyse(Prediction *);
};

extern "C" SensorFrameShift * builder0(int n, DNASeq *X) {  return new SensorFrameShift(n, X); }

#endif
