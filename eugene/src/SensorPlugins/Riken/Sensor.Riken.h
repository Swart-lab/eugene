#ifndef  SENSOR_RIKEN_H_INCLUDED
#define  SENSOR_RIKEN_H_INCLUDED

#include <vector>
#include <algorithm>

#include "../../EuGene/Sensor.h"

// RAFL: Riken Arabidopsis Full Length cDNA
static const int MAX_RIKEN_LENGTH = 60000;
static const int MAX_RIKEN_EST_LENGTH = 3000;
static const int MIN_RIKEN_LENGTH = 120; // 2* riken overlap (60)
static const int MIN_RIKEN_EST_LENGTH = 10;
const REAL RAFLPenalty = -120;

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
class SensorRiken : public Sensor
{
 private:
  std::vector <RAFLgene> RAFL;
  int RAFLpos;                  // Position par rapport a un gene RAFL
  int RAFLindex;                // Index du RIKEN en cours
  int RAFL_A_Traiter;           // Index du RIKEN en cours

 public:
  SensorRiken   (int);
  virtual ~SensorRiken    ();
  virtual void Init       (DNASeq *);
  virtual void GiveInfo   (DNASeq *, int, DATA *);
  virtual void Plot       (DNASeq *);
  virtual void PostAnalyse(Prediction *);
};

extern "C" SensorRiken * builder0( int n ) {  return new SensorRiken(n);}

#endif
