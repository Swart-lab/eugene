#ifndef  SENSOR_NG2_H_INCLUDED
#define  SENSOR_NG2_H_INCLUDED

#include "../../EuGene/Sensor.h"

/*************************************************************
 **                       SensorNG2                         **
 *************************************************************/
class SensorNG2 : public Sensor
{
 private:
  int PositionGiveInfo;

  std::vector<int>  vPosAccF, vPosAccR, vPosDonF, vPosDonR;
  std::vector<REAL> vValAccF, vValAccR, vValDonF, vValDonR;

  int iAccF, iAccR, iDonF, iDonR;
  double accB, accP, donB, donP;
  
  void ReadNG2F(char[FILENAME_MAX+1], int);
  void ReadNG2R(char[FILENAME_MAX+1], int);

 public:
  SensorNG2  (int);
  virtual ~SensorNG2      ();
  virtual void Init       (DNASeq *);
  virtual void GiveInfo   (DNASeq *X, int, DATA *);
  virtual void Plot       (DNASeq *X);
  virtual void PostAnalyse(Prediction *);
};

extern "C" SensorNG2 * builder0( int n ) {  return new SensorNG2(n);}


#endif
