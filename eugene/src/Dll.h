#ifndef __DLL_H
#define __DLL_H

#include "Sensor.h"

class DLLManager
{
 protected:
  void *h;
  const char *err;
  void *(*factory_func)(void);	

 public:
  DLLManager (const char *fname, const char *func_name=0);
  ~DLLManager ();
  bool GetSymbol ( void **, const char *sym_name );
  const char *LastError () { return err; }
};

class DLLFactory : public DLLManager
{
 public:
  SensorFactory *factory;
  
  DLLFactory(const char *fname,
	     const char *func_name=0) : DLLManager( fname, func_name ) {
    if( factory_func )
      factory = (SensorFactory *)factory_func();
    else 
      factory = 0;
  }
  
  ~DLLFactory() { delete factory; }
};

#endif
