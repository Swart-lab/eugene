#ifndef  SENSOR_RIKEN_H_INCLUDED
#define  SENSOR_RIKEN_H_INCLUDED

#include <vector>
#include <algorithm>

#include "../../EuGene/Sensor.h"

// RAFL: Riken Arabidopsis Full Length cDNA
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
  int RAFL_A_Traiter_Remenber;
  int StrandRespect;
  int MIN_EST_DIFF; // default = 100;
  int MAX_OVERLAP; // default = 60;
  int MAX_RIKEN_LENGTH; // default = 60000;
  int MAX_RIKEN_EST_LENGTH; // default  = 3000;
  int MIN_RIKEN_LENGTH; // default  = 120; // 2* riken overlap (60)
  int MIN_RIKEN_EST_LENGTH; // default  = 10;
  REAL RAFLPenalty; // default  = -120;

 public:
  SensorRiken   (int n, DNASeq *X);
  virtual ~SensorRiken    ();
  virtual void Init       (DNASeq *);
  virtual void GiveInfo   (DNASeq *, int, DATA *);
  virtual void Plot       (DNASeq *);
  virtual void PostAnalyse(Prediction *);
};

extern "C" SensorRiken * builder0( int n, DNASeq *X) {  return new SensorRiken(n, X);}

#endif
