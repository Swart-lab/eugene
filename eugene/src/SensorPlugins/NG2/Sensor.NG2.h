#ifndef  SENSOR_NG2_H_INCLUDED
#define  SENSOR_NG2_H_INCLUDED

#include "Sensor.h"

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

#endif
