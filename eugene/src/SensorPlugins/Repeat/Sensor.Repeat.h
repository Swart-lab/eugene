#ifndef  SENSOR_REPEAT_H_INCLUDED
#define  SENSOR_REPEAT_H_INCLUDED

#include "Sensor.h"

/*************************************************************
 **                     SensorRepeat                        **
 *************************************************************/
class SensorRepeat : public Sensor
{
 private:
  unsigned char *ForcedIG;

 public:
  SensorRepeat  ();
  ~SensorRepeat ();
  void Init     (DNASeq *);
  void GiveInfo (DNASeq *, int, DATA *);
};

#endif
