#ifndef __DLL_H
#define __DLL_H

#include "Sensor.h"

class SensorLoader
{
 protected:
  void *h;
  const char *err;
  Sensor *(*builder_func)(int);	

 public:
  SensorLoader (const char *fname, const char *func_name=0);
  ~SensorLoader ();
  Sensor *MakeSensor(int);
  const char *LastError () { return err; }
};

#endif
