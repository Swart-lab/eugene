#ifndef  PARAM_H_INCLUDED
#define  PARAM_H_INCLUDED
#include <cstdio>
#include <map>
#include <strstream>
#include <string>

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

const char VERSION_PAR[FILENAME_MAX+1] = "14_10_02";
const char VERSION[FILENAME_MAX+1]     = "1.46 (141002)";

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
  char *DFT_MATRIX;
  char *DFT_OUTPUT;  
  double FsP;
  char parname[FILENAME_MAX+1];
  int normopt, window, offset, resx, resy, gfrom, gto, golap, glen;
  std::map <const char*, const char*, ltstr> m;
  std::map <const char*, const char*, ltstr>::iterator iter;

  void ReadArg(int, char *[]);
  void ReadPar(char *);
  
 public:
  FILE *fp;
    
  Parameters  ();
  ~Parameters ();
  void   initParam (int, char *[]);
  char*  getC (char *key);
  double getD (char *key);
  int    getI (char *key);
  int    getUseSensor (char **, int*);
  void   set  (const char *key, const char *value);
  const char*  intToChar    (int);
  const char*  doubleToChar (double);
};

#endif
