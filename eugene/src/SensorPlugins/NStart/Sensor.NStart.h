#ifndef  SENSOR_NSTART_H_INCLUDED
#define  SENSOR_NSTART_H_INCLUDED

#include "../../EuGene/Sensor.h"

/*************************************************************
 **                    SensorNStart                         **
 *************************************************************/
class SensorNStart : public Sensor
{
 private:
  std::vector<int>  vPosF, vPosR;
  std::vector<REAL> vValF, vValR;
  std::vector<int>::iterator iter;
  int indexF, indexR;
  double startP, startB;
  
  void ReadNStartF (char[FILENAME_MAX+1], int);
  void ReadNStartR (char[FILENAME_MAX+1], int);

 public:
  SensorNStart   (int);
  virtual ~SensorNStart   ();
  virtual void Init       (DNASeq *);
  virtual void GiveInfo   (DNASeq *X, int, DATA *);
  virtual void Plot       (DNASeq *X);
  virtual void PostAnalyse(Prediction *);
};

extern "C" SensorNStart* builder0( int n ) { return new SensorNStart(n);}

#endif
