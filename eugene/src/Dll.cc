#include <dlfcn.h>
#include "Dll.h"

SensorLoader :: SensorLoader(const char *fname, const char *builder=0)
{
  // Try to open the library now and get any error message.
  h = dlopen( fname, RTLD_NOW );
  err = dlerror();
  
  // Try get the creation function if there is no error yet
  builder_func=0;

  if( LastError()==0 && h) {
    builder_func = (Sensor *(*)())dlsym( h, builder ? builder : "builder0" );
    err=dlerror();
  }
}

SensorLoader :: ~SensorLoader()
{
  // Close the library if it isn't null
  if( h != 0 )
    dlclose(h);
}

Sensor* SensorLoader :: MakeSensor()
{
  return (* builder_func)();
}
