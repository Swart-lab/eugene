#ifndef  SENSOR_MARKOVIMM_H_INCLUDED
#define  SENSOR_MARKOVIMM_H_INCLUDED

#include "../../EuGene/Sensor.h"
#include "../0_SensorTk/BStrArray.h"

/****************************************************************
 **                     SensorMarkovIMM                        **
 ****************************************************************/
class SensorMarkovIMM : public Sensor
{
 private:
  BString_Array *IMMatrix[7];
  double transCodant, transIntron, transInter, transUTR5, transUTR3;
  double minGC,maxGC;
  const static int  MODEL_LEN = 9;
  const static int  SIMPLE_MODEL_LEN = 6;
  const static int  ALPHABET_SIZE = 4;

 public:
  SensorMarkovIMM  (int);
  virtual ~SensorMarkovIMM   ();
  virtual void Init       (DNASeq *);
  virtual void ResetIter  ();
  virtual void GiveInfo   (DNASeq *, int, DATA *);
  virtual void GiveInfoAt (DNASeq *, int, DATA *);
  virtual void Plot(DNASeq *);
  virtual void PostAnalyse(Prediction *);
};

extern "C" SensorMarkovIMM * builder0( int n ) { return new SensorMarkovIMM(n);}

#endif
