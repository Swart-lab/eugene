#ifndef  SENSOR_EST_H_INCLUDED
#define  SENSOR_EST_H_INCLUDED

#include "../Sensor.h"
#include "../Hits.h"

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
  double estP, estM;
  
  Hits** ESTAnalyzer (FILE *, unsigned char *, int, int *, DNASeq *);
  void   ESTSupport  (Prediction *pred, int Tdebut, int Tfin,
		      int debut,int fin,  Hits **HitTable, int Size);

 public:
  SensorEst  (int);
  virtual ~SensorEst      ();
  virtual void Init       (DNASeq *);
  virtual void ResetIter  ();
  virtual void GiveInfo   (DNASeq *, int, DATA *);
  virtual void GiveInfoAt (DNASeq *, int, DATA *);
  virtual void Plot       (DNASeq *);
  virtual void PostAnalyse(Prediction *);
};

extern "C" SensorEst * builder0( int n ) { return new SensorEst(n);}

#endif
