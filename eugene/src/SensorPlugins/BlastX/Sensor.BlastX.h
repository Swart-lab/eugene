#ifndef  SENSOR_BLASTX_H_INCLUDED
#define  SENSOR_BLASTX_H_INCLUDED

#include "../Sensor.h"

/*************************************************************
 **                     SensorBlastX                        **
 *************************************************************/
class SensorBlastX : public Sensor
{
 private:
  REAL *ProtMatch, *ProtMatchLevel;
  int  *ProtMatchPhase;
  double keyBXLevel[8];
  int    minL8;
  
  void LoadContentScore (DNASeq *);
  int  PhaseAdapt (char);
  char ph06       (char);

 public:
  SensorBlastX  (int);
  ~SensorBlastX ();
  void Init     (DNASeq *);
  void GiveInfo (DNASeq *, int, DATA *);
};

extern "C" SensorBlastX * builder0( int n ) {  return new SensorBlastX(n);}
#endif
