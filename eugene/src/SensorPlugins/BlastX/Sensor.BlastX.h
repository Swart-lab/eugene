#ifndef  SENSOR_BLASTX_H_INCLUDED
#define  SENSOR_BLASTX_H_INCLUDED

#include "../Sensor.h"

/*************************************************************
 **                     SensorBlastX                        **
 *************************************************************/
class SensorBlastX : public Sensor
{
 private:
  REAL *ProtMatch, *ProtMatchLevel;
  int  *ProtMatchPhase;
  double keyBXLevel[8];
  int    minL8;
  
  void LoadContentScore (DNASeq *);
  int  PhaseAdapt (char);
  char ph06       (char);

 public:
  SensorBlastX  ();
  ~SensorBlastX ();
  void Init     (DNASeq *);
  void GiveInfo (DNASeq *, int, DATA *);
};

/*************************************************************
 **                  SensorBlastXFactory                    **
 *************************************************************/
class SensorBlastXFactory : public SensorFactory
{
 public:
  SensorBlastXFactory()  { }
  
  ~SensorBlastXFactory() { }
  
  virtual Sensor * CreateSensor() {
    return new SensorBlastX;
  }
};

extern "C" void * factory0( void ) {
  return new SensorBlastXFactory;
}

#endif
