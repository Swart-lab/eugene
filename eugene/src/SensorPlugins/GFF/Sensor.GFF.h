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
  char  name[30];
  char  feature[30];
  int   start, end;
  char  strand, frame;
  int   A, T, C, G, N;
  float GC;
  
  GFFObject  ();
  GFFObject  (char*, char*,int,int,char,char,int,int,int,int,int,float);
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
  static FILE *ppfile;
  static bool IsInitialized;

  void ReadGFF(char[FILENAME_MAX+1]);

 public:
  SensorGFF          (int n, DNASeq *X);
  virtual ~SensorGFF ();
  virtual void Init       (DNASeq *);
  virtual void GiveInfo   (DNASeq *X, int, DATA *);
  virtual void Plot       (DNASeq *X);
  virtual void PostAnalyse(Prediction *);
};

extern "C" SensorGFF* builder0( int n, DNASeq *X) { return new SensorGFF(n, X);}



FILE * SensorGFF::ppfile;
bool SensorGFF::IsInitialized = false;


#endif
