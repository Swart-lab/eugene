#ifndef  PARAM_H_INCLUDED
#define  PARAM_H_INCLUDED

#include <map>

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <string.h>
#ifdef HAVE_STRINGS_H
#include <strings.h>
#endif

#ifdef HAVE_GETOPT_H
#include <getopt.h>
#else
#include "getopt.h"
#endif
#include "Const.h"
#include "System.h"

const char VERSION_PAR[FILENAME_MAX+1] = "28_08_03";
const char VERSION[FILENAME_MAX+1]     = "1.58 (280803)";

extern char   *optarg;   
extern int     optind;


class ltstr
{
  public :
    bool operator() (const char* s1, const char* s2) const
    { return (strcmp(s1, s2) <0); };
};

class Parameters
{
 private: 
  std::map <const char*, const char*, ltstr> m;
  std::map <const char*, const char*, ltstr>::iterator iter;

  void ReadArg(int, char *[]);
  void ReadPar(char *);
  
 public:
  FILE *fp;
    
  Parameters  ();
  ~Parameters ();
  void   initParam (int, char *[]);
  int    count(char *key);
  char*  getC (char *key, int index = 0);
  double getD (char *key, int index = 0);
  int    getI (char *key, int index = 0);
  int    getUseSensor (char **, int*);
  void   set  (const char *key, const char *value);
};

#endif
