#include <dlfcn.h>
#include "Dll.h"

DLLManager :: DLLManager(const char *fname, const char *factory=0)
{
  // Try to open the library now and get any error message.
  h = dlopen( fname, RTLD_NOW );
  err = dlerror();
  
  // Try get the factory function if there is no error yet
  factory_func=0;

  if( LastError()==0 ) {		
    GetSymbol( (void **)&factory_func, factory ? factory : "factory0" );
  } 
}

DLLManager :: ~DLLManager()
{
  // Close the library if it isn't null
  if( h != 0 )
    dlclose(h);
}

bool DLLManager :: GetSymbol(void **v, const char *sym_name)
{
  // Try extract a symbol from the library get any error message
  // is there is any.
  if( h!=0 ) {
    *v = dlsym( h, sym_name );
    err=dlerror();
    if( err==0 )
      return true;
    else
      return false;
  }
  else
    return false;
}
