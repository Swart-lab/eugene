#ifndef  SENSOR_EST_H_INCLUDED
#define  SENSOR_EST_H_INCLUDED

#include "Sensor.h"
#include "Hits.h"

/*************************************************************
 **                       SensorEst                         **
 *************************************************************/
class SensorEst : public Sensor
{
 private:
  unsigned char *ESTMatch;
  Hits **HitTable;
  int NumEST;

 public:
  SensorEst  ();
  ~SensorEst ();
  void Init     (DNASeq *);
  void GiveInfo (DNASeq *, int, DATA *);
  Hits** ESTAnalyzer (FILE *, unsigned char *, int, int *, DNASeq *);
};

#endif
