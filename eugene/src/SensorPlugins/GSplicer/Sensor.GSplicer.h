#ifndef  SENSOR_GSPLICER_H_INCLUDED
#define  SENSOR_GSPLICER_H_INCLUDED

#include "../../EuGene/Sensor.h"

/*************************************************************
 **                      SensorGSplicer
 *************************************************************/
class SensorGSplicer : public Sensor
{
 private:
  std::vector<int>  vPosAccF, vPosAccR, vPosDonF, vPosDonR;
  std::vector<REAL> vValAccF, vValAccR, vValDonF, vValDonR;
  std::vector<int>::iterator iter;
  int iAccF, iAccR, iDonF, iDonR;
  double coefAcc, penAcc, coefDon, penDon;
  
  void ReadGSplicer(char[FILENAME_MAX+1], int);
  
 public:
  SensorGSplicer          (int);
  virtual ~SensorGSplicer ();
  virtual void Init       (DNASeq *);
  virtual void GiveInfo   (DNASeq *X, int, DATA *);
  virtual void Plot       (DNASeq *X);
  virtual void PostAnalyse(Prediction *);
};

extern "C" SensorGSplicer* builder0( int n ) { return new SensorGSplicer(n);}

#endif
