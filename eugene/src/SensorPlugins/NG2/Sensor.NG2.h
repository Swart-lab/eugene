#ifndef  SENSOR_NG2_H_INCLUDED
#define  SENSOR_NG2_H_INCLUDED

#include "../Sensor.h"

/*************************************************************
 **                       SensorNG2                         **
 *************************************************************/
class SensorNG2 : public Sensor
{
 private:
  REAL *Acc[2];
  REAL *Don[2];
  double accB, accP, donB, donP;

  void ReadNG2(char[FILENAME_MAX+1], int, int);

 public:
  SensorNG2  ();
  ~SensorNG2 ();
  void Init     (DNASeq *);
  void GiveInfo (DNASeq *, int, DATA *);
};

/*************************************************************
 **                   SensorNG2Factory                      **
 *************************************************************/
class SensorNG2Factory : public SensorFactory
{
public:
  SensorNG2Factory()  { }
	
  ~SensorNG2Factory() { }
  
  virtual Sensor * CreateSensor() {
    return new SensorNG2;
  }
};

extern "C" void * factory0( void ) {
  return new SensorNG2Factory;
}

#endif
