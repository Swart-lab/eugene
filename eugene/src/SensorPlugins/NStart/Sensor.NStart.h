#ifndef  SENSOR_NSTART_H_INCLUDED
#define  SENSOR_NSTART_H_INCLUDED

#include "../Sensor.h"

/*************************************************************
 **                    SensorNStart                         **
 *************************************************************/
class SensorNStart : public Sensor
{
 private:
  std::vector<int>  vPosF, vPosR;
  std::vector<REAL> vValF, vValR;
  std::vector<int>::iterator iter;
  int iterF, iterR;
  double startP, startB;
  
  void ReadNStartF (char[FILENAME_MAX+1], int);
  void ReadNStartR (char[FILENAME_MAX+1], int);

 public:
  SensorNStart   (int);
  virtual ~SensorNStart   ();
  virtual void Init       (DNASeq *);
  virtual void ResetIter  ();
  virtual void GiveInfo   (DNASeq *, int, DATA *);
  virtual void GiveInfoAt (DNASeq *, int, DATA *);
};

extern "C" SensorNStart* builder0( int n ) { return new SensorNStart(n);}

#endif
