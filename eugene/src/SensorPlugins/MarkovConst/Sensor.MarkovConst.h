#ifndef  SENSOR_MARKOVCONST_INCLUDED
#define  SENSOR_MARKOVCONST_INCLUDED

#include "../Sensor.h"

/*************************************************************
 **                     SensorMarkovConst                    **
 *************************************************************/
class SensorMarkovConst : public Sensor
{
 private:
  REAL value;
  double transCodant, transIntron, transInter, transUTR5, transUTR3;
  double minGC,maxGC;

 public:
  SensorMarkovConst  (int);
  virtual ~SensorMarkovConst   ();
  virtual void Init       (DNASeq *);
  virtual void GiveInfo   (DNASeq *, int, DATA *);
  virtual void Plot       (DNASeq *);
  virtual void PostAnalyse(Prediction *);
};

extern "C" SensorMarkovConst * builder0(int n) {  return new SensorMarkovConst(n); }

#endif
