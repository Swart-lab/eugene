#ifndef  SENSOR_SPred_H_INCLUDED
#define  SENSOR_SPred_H_INCLUDED

#include "../Sensor.h"

/*************************************************************
 **                     SensorSPred                         **
 *************************************************************/
class SensorSPred : public Sensor
{
 private:
  std::vector<int>  vPosAccF, vPosAccR, vPosDonF, vPosDonR;
  std::vector<REAL> vValAccF, vValAccR, vValDonF, vValDonR;
  std::vector<int>::iterator iter;
  int iAccF, iAccR, iDonF, iDonR;
  double accP, accB, donP, donB;
  
  void ReadSPredF(char[FILENAME_MAX+1], int);
  void ReadSPredR(char[FILENAME_MAX+1], int);
  
 public:
  SensorSPred   (int);
  virtual ~SensorSPred    ();
  virtual void Init       (DNASeq *);
  virtual void GiveInfo   (DNASeq *X, int, DATA *);
  virtual void Plot       (DNASeq *X);
  virtual void PostAnalyse(Prediction *);
};

extern "C" SensorSPred * builder0( int n ) {  return new SensorSPred(n);}

#endif
