#ifndef  SENSOR_RIKEN_H_INCLUDED
#define  SENSOR_RIKEN_H_INCLUDED

#include <vector>
#include <algorithm>

#include "../Sensor.h"

/*************************************************************
 **                     SensorRiken                         **
 *************************************************************/
// ************
// * RAFLgene *     // RAFL: Riken Arabidopsis Full Length cDNA 
// ************
class RAFLgene
{
 private:

 public:
  int deb;
  int fin;
  signed char sens;
  char ID[FILENAME_MAX+1];

  RAFLgene  ();
  ~RAFLgene ();
};

// ****************
// * Sensor Riken *
// ****************
const REAL RAFLPenalty = NINFINITY; 

class SensorRiken : public Sensor
{
 private:
  std::vector <RAFLgene*> RAFL;
  int RAFLpos;                  // Position par rapport a un gene RAFL
  int RAFLindex;                // Index du RIKEN en cours
  int RAFL_A_Traiter;           // Index du RIKEN en cours

 public:
  SensorRiken   (int);
  virtual ~SensorRiken    ();
  virtual void Init       (DNASeq *);
  virtual void ResetIter  ();
  virtual void GiveInfo   (DNASeq *, int, DATA *);
  virtual void GiveInfoAt (DNASeq *, int, DATA *);
  virtual void Plot       (DNASeq *);
  virtual void PostAnalyse(Prediction *);
};

extern "C" SensorRiken * builder0( int n ) {  return new SensorRiken(n);}

#endif
