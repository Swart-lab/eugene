#ifndef  SENSOR_REPEAT_H_INCLUDED
#define  SENSOR_REPEAT_H_INCLUDED

#include "../Sensor.h"

/*************************************************************
 **                     SensorRepeat                        **
 *************************************************************/
class SensorRepeat : public Sensor
{
 private:
  std::vector<int>           vPos;
  std::vector<int>::iterator iter;
  int index;

 public:
  SensorRepeat  (int);
  virtual ~SensorRepeat   ();
  virtual void Init       (DNASeq *);
  virtual void ResetIter  ();
  virtual void GiveInfo   (DNASeq *, int, DATA *);
  virtual void GiveInfoAt (DNASeq *, int, DATA *);
};

extern "C" SensorRepeat * builder0( int n ) { return new SensorRepeat(n);}

#endif
