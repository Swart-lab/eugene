#ifndef  SENSOR_SPred_H_INCLUDED
#define  SENSOR_SPred_H_INCLUDED

#include "../Sensor.h"

/*************************************************************
 **                     SensorSPred                         **
 *************************************************************/
class SensorSPred : public Sensor
{
 private:
  REAL *Acc[2];
  REAL *Don[2];
  double accP, accB, donP, donB;

  void ReadSPred(char[FILENAME_MAX+1], int, int);

 public:
  SensorSPred   ();
  ~SensorSPred  ();
  void Init     (DNASeq *);
  void GiveInfo (DNASeq *, int, DATA *);
};

extern "C" SensorSPred * builder0( void ) {  return new SensorSPred;}

#endif
