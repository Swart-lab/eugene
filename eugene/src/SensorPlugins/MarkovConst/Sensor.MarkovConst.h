#ifndef  SENSOR_MARKOVCONST_INCLUDED
#define  SENSOR_MARKOVCONST_INCLUDED

#include "../../EuGene/Sensor.h"

/*************************************************************
 **                     SensorMarkovConst                    **
 *************************************************************/
class SensorMarkovConst : public Sensor
{
 private:
  double transCodant, transIntron, transIntronU, transInter, transUTR5, transUTR3;
  double minGC,maxGC;

 public:
  SensorMarkovConst  (int n, DNASeq *X);
  virtual ~SensorMarkovConst   ();
  virtual void Init       (DNASeq *);
  virtual void GiveInfo   (DNASeq *, int, DATA *);
  virtual void Plot       (DNASeq *);
  virtual void PostAnalyse(Prediction *);
};

extern "C" SensorMarkovConst * builder0(int n, DNASeq *X) {  return new SensorMarkovConst(n, X); }

#endif
