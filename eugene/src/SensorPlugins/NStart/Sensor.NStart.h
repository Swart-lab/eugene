#ifndef  SENSOR_NSTART_H_INCLUDED
#define  SENSOR_NSTART_H_INCLUDED

#include "../Sensor.h"

/*************************************************************
 **                    SensorNStart                         **
 *************************************************************/
class SensorNStart : public Sensor
{
 private:
  REAL *Start[2];
  double startP, startB;
  
  void ReadNStart (char[FILENAME_MAX+1], int, int);

 public:
  SensorNStart   ();
  ~SensorNStart  ();
  void Init     (DNASeq *);
  void GiveInfo (DNASeq *, int, DATA *);
};


/*************************************************************
 **                 SensorNStartFactory                     **
 *************************************************************/
class SensorNStartFactory : public SensorFactory
{
 public:
  SensorNStartFactory()  { }
  
  ~SensorNStartFactory() { }
  
  virtual Sensor * CreateSensor() {
    return new SensorNStart;
  }
};

extern "C" void * factory0( void ) {
  return new SensorNStartFactory;
}

#endif
