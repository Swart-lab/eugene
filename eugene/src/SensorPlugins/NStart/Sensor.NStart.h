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
  SensorNStart   (int);
  ~SensorNStart  ();
  void Init     (DNASeq *);
  void GiveInfo (DNASeq *, int, DATA *);
};

extern "C" SensorNStart* builder0( int n ) { return new SensorNStart(n);}

#endif
