#ifndef  SENSOR_EST_H_INCLUDED
#define  SENSOR_EST_H_INCLUDED

#include "../Sensor.h"
#include "../Hits.h"

/*************************************************************
 **                       SensorEst                         **
 *************************************************************/
class SensorEst : public Sensor
{
 private:
  unsigned char *ESTMatch;
  Hits **HitTable;
  int NumEST;
  double estP, estM;

 public:
  SensorEst  ();
  ~SensorEst ();
  void Init     (DNASeq *);
  void GiveInfo (DNASeq *, int, DATA *);
  Hits** ESTAnalyzer (FILE *, unsigned char *, int, int *, DNASeq *);
};

/*************************************************************
 **                    SensorEstFactory                     **
 *************************************************************/
class SensorEstFactory : public SensorFactory
{
 public:
  SensorEstFactory()  { }
  
  ~SensorEstFactory() { }
  
  virtual Sensor * CreateSensor() {
    return new SensorEst;
  }
};

extern "C" void * factory0( void ) {
  return new SensorEstFactory;
}

#endif