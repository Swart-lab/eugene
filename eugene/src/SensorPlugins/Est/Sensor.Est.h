#ifndef  SENSOR_EST_H_INCLUDED
#define  SENSOR_EST_H_INCLUDED

#include "../../EuGene/Sensor.h"
#include "../../EuGene/Hits.h"

/*************************************************************
 **                       SensorEst                         **
 *************************************************************/
class SensorEst : public Sensor
{
 private:
  std::vector<int>  vPos;
  std::vector<unsigned char> vESTMatch;
  std::vector<int>::iterator iter;
  int index;
  unsigned char *ESTMatch;
  Hits **HitTable;
  int NumEST;
  double estP, utrP;
  double DonorThreshold;
  int    estM, utrM;
  int N;
  
  Hits** ESTAnalyzer (FILE *, unsigned char *, int, int *, DNASeq *);
  void   ESTSupport  (Prediction *pred, int Tdebut, int Tfin,
		      int debut,int fin,  Hits **HitTable, int Size);

 public:
  SensorEst  (int n, DNASeq *X);
  virtual ~SensorEst      ();
  virtual void Init       (DNASeq *);
  virtual void GiveInfo   (DNASeq *, int, DATA *);
  virtual void Plot       (DNASeq *);
  virtual void PostAnalyse(Prediction *);
};

extern "C" SensorEst * builder0( int n, DNASeq *X) { return new SensorEst(n, X);}

#endif
