#ifndef  SENSOR_BLASTX_H_INCLUDED
#define  SENSOR_BLASTX_H_INCLUDED

#include "../../EuGene/Sensor.h"

/*************************************************************
 **                     SensorBlastX                        **
 *************************************************************/
class SensorBlastX : public Sensor
{
 private:
  REAL *ProtMatch, *ProtMatchLevel;
  int  *ProtMatchPhase;
  std::vector<int>  vPos,     vPMPhase;
  std::vector<REAL> vPMLevel, vPMatch;
  std::vector<int>::iterator iter;
  int    index;
  double keyBXLevel[8];
  int    minL8;
  
  void LoadContentScore (DNASeq *);
  int  PhaseAdapt (char);
  char ph06       (char);

 public:
  SensorBlastX  (int);
  virtual ~SensorBlastX   ();
  virtual void Init       (DNASeq *);
  virtual void GiveInfo   (DNASeq *, int, DATA *);
  virtual void Plot       (DNASeq *);
  virtual void PostAnalyse(Prediction *);
};

extern "C" SensorBlastX * builder0( int n ) {  return new SensorBlastX(n);}
#endif
