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
  static BString_Array *IMMatrix[7];
  static double minGC;
  static double maxGC;
  static bool IsInitialized;
  const static int  MODEL_LEN = 9;
  const static int  SIMPLE_MODEL_LEN = 6;
  const static int  ALPHABET_SIZE = 4;

 public:
  SensorMarkovIMM  (int n, DNASeq *X);
  virtual ~SensorMarkovIMM   ();
  virtual void Init       (DNASeq *);
  virtual void ResetIter  ();
  virtual void GiveInfo   (DNASeq *, int, DATA *);
  virtual void GiveInfoAt (DNASeq *, int, DATA *);
  virtual void Plot(DNASeq *);
  virtual void PostAnalyse(Prediction *);
};

extern "C" SensorMarkovIMM * builder0( int n, DNASeq *X) { return new SensorMarkovIMM(n, X);}



BString_Array* SensorMarkovIMM::IMMatrix[7];
double SensorMarkovIMM::minGC;
double SensorMarkovIMM::maxGC;
bool SensorMarkovIMM::IsInitialized = false;


#endif
