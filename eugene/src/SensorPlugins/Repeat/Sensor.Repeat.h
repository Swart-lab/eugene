#ifndef  SENSOR_REPEAT_H_INCLUDED
#define  SENSOR_REPEAT_H_INCLUDED

#include "../Sensor.h"

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

/*************************************************************
 **                   SensorRepeatFactory                   **
 *************************************************************/
class SensorRepeatFactory : public SensorFactory
{
 public:
  SensorRepeatFactory()  { }
	
  ~SensorRepeatFactory() { }
  
  virtual Sensor * CreateSensor() {
    return new SensorRepeat;
  }
};

extern "C" void * factory0( void ) {
  return new SensorRepeatFactory;
}

#endif
