#ifndef  SENSOR_GCPLOT_H_INCLUDED
#define  SENSOR_GCPLOT_H_INCLUDED

#include "../../EuGene/Sensor.h"

/*************************************************************
 **                      SensorGCPlot
 *************************************************************/
class SensorGCPlot : public Sensor
{
 private:
  int Color;
  int Window;
  int Up;
  int Over;
  double Zoom;
  double Zoom3;

 public:
  SensorGCPlot          (int n, DNASeq *X);
  virtual ~SensorGCPlot ();
  virtual void Init       (DNASeq *X);
  virtual void GiveInfo   (DNASeq *X, int, DATA *);
  virtual void Plot       (DNASeq *X);
  virtual void PostAnalyse(Prediction *);
};

extern "C" SensorGCPlot* builder0( int n, DNASeq *X) { return new SensorGCPlot(n, X);}

#endif
