#ifndef  PARAM_H_INCLUDED
#define  PARAM_H_INCLUDED

#include <map>
#include <vector>

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <string>
#include <string.h>
#ifdef HAVE_STRINGS_H
#include <strings.h>
#endif

#ifdef HAVE_GETOPT_H
#include <getopt.h>
#else
#ifndef HAVE_GETOPT
#include "getopt.h"
#endif
#endif



#include "Const.h"
#include "System.h"

const char VERSION[FILENAME_MAX+1]     = "2.0";
const char VERSION_DATE[FILENAME_MAX+1]     = "8080704";
const char VERSION_PAR[FILENAME_MAX+1]     = "8080704";
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
    
  ~Parameters ();
  void   initParam (int, char *[]);
  int    count(char *key);
  char*  getC (char *key, int index = 0);
  double getD (char *key, int index = 0);
  int    getI (char *key, int index = 0);
  int    getUseSensor (char **, int*);
  void   set  (const char *key, const char *value);
  void   setD (const char *key, double n);
  std::string WriteParam (const char* para_file, std::vector<std::string> para_name, 
			  std::vector<double> para_val);

  void ResetIter(void);
};

#endif
