#ifndef  MSENSOR_H_INCLUDED
#define  MSENSOR_H_INCLUDED
#include <cstdio>
#include <vector>
#include <algorithm>

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <string.h>
#ifdef HAVE_STRINGS_H
#include <strings.h>
#endif

#include "Sensor.h"
#include "DNASeq.h"
#include "Param.h"
#include "Const.h"
#include "Dll.h"

/*************************************************************
 **                        UseSensor                        **
 *************************************************************/
class UseSensor
{
 private:
  
 public:
  int  Priority;
  char Name[FILENAME_MAX+1];

  UseSensor  ();
  UseSensor  (int, char[FILENAME_MAX+1]);
  ~UseSensor ();
};

/*************************************************************
 **                      MasterSensor                       **
 *************************************************************/
class MasterSensor
{
 private:
  std::vector <UseSensor*> msList;
  std::vector <Sensor*>    theSensors;
  char **soList;
  char **useList;
  SensorLoader **dllList;
  
 public:
  ~MasterSensor (void);
  void InitMaster   (DNASeq *X);
  void InitSensors  (DNASeq *X);
  void GetInfoAt    (DNASeq *X, int pos, DATA *d);
  void PrintDataAt  (DNASeq *X, int pos, DATA *d);
  int  GetInfoSpAt  (TYPE_SENSOR type, DNASeq *X, int pos, DATA *d);
  void PostAnalyse  (Prediction *pred);
};

#endif
