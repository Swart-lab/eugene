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
  MasterSensor  ();
  ~MasterSensor ();
  void InitMaster   ();
  void InitSensors  (DNASeq *);
  void ResetIterator();
  void ResetType    ();
  void GetInfoAt    (DNASeq *, int, DATA *);
  int  GetInfoSpAt  (TYPE_SENSOR, DNASeq *, int, DATA *);
  void ResetSensors ();
};

#endif
