#ifndef  SENSOR_BLASTX_H_INCLUDED
#define  SENSOR_BLASTX_H_INCLUDED

#include "../../EuGene/Sensor.h"

/*************************************************************
 **                     SensorBlastX                        **
 *************************************************************/
class SensorBlastX : public Sensor
{
 private:
  double *ProtMatch, *ProtMatchLevel;
  int  *ProtMatchPhase;
  std::vector<int>  vPos,     vPMPhase;
  std::vector<double> vPMLevel, vPMatch;
  std::vector<int>::iterator iter;
  int    index;
  double keyBXLevel[10];
  int    minIn;
  int blastxM;
  int N;

  void LoadContentScore (DNASeq *);
  int  PhaseAdapt (char);
  char ph06       (char);

 public:
  SensorBlastX  (int n, DNASeq *X);
  virtual ~SensorBlastX   ();
  virtual void Init       (DNASeq *);
  virtual void GiveInfo   (DNASeq *, int, DATA *);
  virtual void Plot       (DNASeq *);
  virtual void PostAnalyse(Prediction *);
};

extern "C" SensorBlastX * builder0( int n, DNASeq *X) {  return new SensorBlastX(n, X);}
#endif
