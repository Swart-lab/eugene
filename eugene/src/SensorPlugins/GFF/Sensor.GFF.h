#ifndef  SENSOR_GFF_H_INCLUDED
#define  SENSOR_GFF_H_INCLUDED

#include "../../EuGene/Sensor.h"

/*************************************************************
 **                        GFFObject                        **
 *************************************************************/
class GFFObject
{
 private:
  
 public:
  char feature[30];
  int  start, end;
  char strand, frame;

  GFFObject  ();
  GFFObject  (char*,int,int,char,char);
  void Print       ();
  void PrintHeader ();
  ~GFFObject ();
};

/*************************************************************
 **                      SensorGFF
 *************************************************************/
class SensorGFF : public Sensor
{
 private:
  std::vector <GFFObject*> gffList;
  void ReadGFF(char[FILENAME_MAX+1]);
  
 public:
  SensorGFF          (int);
  virtual ~SensorGFF ();
  virtual void Init       (DNASeq *);
  virtual void GiveInfo   (DNASeq *X, int, DATA *);
  virtual void Plot       (DNASeq *X);
  virtual void PostAnalyse(Prediction *);
};

extern "C" SensorGFF* builder0( int n ) { return new SensorGFF(n);}

#endif
