#ifndef  SENSOR_STARTCONST_INCLUDED
#define  SENSOR_STARTCONST_INCLUDED

#include "../../EuGene/Sensor.h"

/*************************************************************
 **                     SensorTranscript                    **
 *************************************************************/
class SensorTranscript : public Sensor
{
 private:
  // proba. of transcription Start/Stop
  double transStart;
  double transStop;

 public:
  SensorTranscript  (int n, DNASeq *X);
  virtual ~SensorTranscript   ();
  virtual void Init       (DNASeq *);
  virtual void GiveInfo   (DNASeq *X, int, DATA *);
  virtual void Plot       (DNASeq *X);
  virtual void PostAnalyse(Prediction *);
};

extern "C" SensorTranscript * builder0(int n, DNASeq *X) {  return new SensorTranscript(n, X); }

#endif
