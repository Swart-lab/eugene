#ifndef  SENSOR_MARKOVPROT_INCLUDED
#define  SENSOR_MARKOVPROT_INCLUDED

#include "../../EuGene/Sensor.h"
#include "../0_SensorTk/BStrArray.h"
#include "../0_SensorTk/markov.h"
#include "../0_SensorTk/markov.cc"

/*************************************************************
 **                     SensorMarkovProt                        **
 *************************************************************/
class SensorMarkovProt : public Sensor
{
 private:
  int maxorder, order;
//  char[4] type;
//  TabChaine<ChainePROT21,double>* ModeleProt;
  TabChaine<ChainePROT21,unsigned short>* ModeleProt; 
  TabChaine<ChaineADN,double>* Probacodon;
  double GCrate;
  double transCodant, transIntron, transInter, transUTR5, transUTR3;
  double minGC,maxGC;

 public:
  SensorMarkovProt  (int);
  virtual ~SensorMarkovProt   ();
  virtual void Init       (DNASeq *);
  virtual void ResetIter  ();
  virtual void GiveInfo   (DNASeq *, int, DATA *);
  virtual void GiveInfoAt (DNASeq *, int, DATA *);
  virtual void Plot(DNASeq *);
  virtual void PostAnalyse(Prediction *);
};

extern "C" SensorMarkovProt * builder0( int n ) { return new SensorMarkovProt(n);}

#endif
