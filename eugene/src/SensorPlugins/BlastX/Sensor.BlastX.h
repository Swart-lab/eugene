#ifndef  SENSOR_BLASTX_H_INCLUDED
#define  SENSOR_BLASTX_H_INCLUDED

#include "Sensor.h"

/*************************************************************
 **                     SensorBlastX                        **
 *************************************************************/
class SensorBlastX : public Sensor
{
 private:
  REAL *ProtMatch, *ProtMatchLevel;
  int  *ProtMatchPhase;
  void LoadContentScore (DNASeq *);
  int  PhaseAdapt (char);
  char ph06       (char);

 public:
  SensorBlastX  ();
  ~SensorBlastX ();
  void Init     (DNASeq *);
  void GiveInfo (DNASeq *, int, DATA *);
};

#endif
