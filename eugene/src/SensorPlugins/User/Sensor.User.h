#ifndef  SENSOR_USER_H_INCLUDED
#define  SENSOR_USER_H_INCLUDED

#include "../Sensor.h"
#include "../structure.h"
#include "../yacc.tab.h"

/*************************************************************
 **                       SensorUser                        **
 *************************************************************/
class SensorUser : public Sensor
{
 private:

 public:
  SensorUser  (int);
  ~SensorUser ();
  void Init     (DNASeq *);
  void GiveInfo (DNASeq *, int, DATA *);
};

extern "C" SensorUser * builder0( int n ) {  return new SensorUser(n);}

#endif
