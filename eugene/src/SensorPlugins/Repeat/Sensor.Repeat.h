#ifndef  SENSOR_REPEAT_H_INCLUDED
#define  SENSOR_REPEAT_H_INCLUDED

#include "../../EuGene/Sensor.h"


/*************************************************************
 **                     SensorRepeat                        **
 *************************************************************/
class SensorRepeat : public Sensor
{
 private:
  int PositionGiveInfo;

  std::vector<int>           vDeb;
  std::vector<int>           vFin;
  int index;
  double intronPenalty;
  double UTRPenalty;
  double exonPenalty;

 public:
  SensorRepeat  (int);
  virtual ~SensorRepeat   ();
  virtual void Init       (DNASeq *);
  virtual void GiveInfo   (DNASeq *, int, DATA *);
  virtual void Plot       (DNASeq *);
  virtual void PostAnalyse(Prediction *);
};

extern "C" SensorRepeat * builder0( int n ) { return new SensorRepeat(n);}

#endif
