#include <dlfcn.h>
#include "Dll.h"

SensorLoader :: SensorLoader(const char *fname, const char *builder)
{
  // Try to open the library now and get any error message.
  h = dlopen( fname, RTLD_NOW );
  err = dlerror();
  if (err) printf("dlerror: %s\n",err);
  
  // Try get the creation function if there is no error yet
  builder_func=0;

  if( LastError()==0 && h) {
    builder_func = (Sensor *(*)(int, DNASeq *))dlsym( h, builder ? builder : "builder0" );
    err=dlerror();
  }
}

SensorLoader :: ~SensorLoader()
{
  // Close the library if it isn't null
  if( h != 0 )
    dlclose(h);
}

Sensor* SensorLoader :: MakeSensor(int n, DNASeq *X)
{
  return (* builder_func)(n, X);
}
