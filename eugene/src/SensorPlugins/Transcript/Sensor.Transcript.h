#ifndef  SENSOR_STARTCONST_INCLUDED
#define  SENSOR_STARTCONST_INCLUDED

#include "../Sensor.h"

/*************************************************************
 **                     SensorTranscript                    **
 *************************************************************/
class SensorTranscript : public Sensor
{
 private:
  // proba. of transcription Start/Stop
  REAL transStart;
  REAL transStop;

 public:
  SensorTranscript  (int);
  virtual ~SensorTranscript   ();
  virtual void Init       (DNASeq *);
  virtual void GiveInfo   (DNASeq *X, int, DATA *);
  virtual void Plot       (DNASeq *X);
  virtual void PostAnalyse(Prediction *);
};

extern "C" SensorTranscript * builder0(int n) {  return new SensorTranscript(n); }

#endif
