#ifndef  SENSOR_EU_H_INCLUDED
#define  SENSOR_EU_H_INCLUDED

#include "../Sensor.h"

/*************************************************************
 **                     SensorEuStop                        **
 *************************************************************/
class SensorEuStop : public Sensor
{
 private:
  std::vector<int>  vPosF, vPosR;
  std::vector<REAL> vValF, vValR;
  std::vector<int>::iterator iter;
  int iterF, iterR;

 public:
  SensorEuStop  (int);
  virtual ~SensorEuStop   ();
  virtual void Init       (DNASeq *);
  virtual void ResetIter  ();
  virtual void GiveInfo   (DNASeq *, int, DATA *);
  virtual void GiveInfoAt (DNASeq *, int, DATA *);
};

extern "C" SensorEuStop * builder0(int n) {  return new SensorEuStop(n); }

#endif
