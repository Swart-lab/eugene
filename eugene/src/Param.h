#ifndef  PARAM_H_INCLUDED
#define  PARAM_H_INCLUDED
#include <cstdio>

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

const char VERSION_PAR[FILENAME_MAX+1] = "21_08_02";
const char VERSION[FILENAME_MAX+1]     = "1.2b (280802)";

extern char   *optarg;   
extern int     optind;

class Parameters
{
 private:
  char *DFT_MATRIX;
  char *DFT_OUTPUT;

  void ReadArg(int, char *[]);
  void ReadPar(char *);

 public:
  // Arguments
  REAL ExonPrior, IntronPrior, InterPrior,FivePrimePrior, ThreePrimePrior;
  int normopt, blastopt, estopt, estanal, ncopt, raflopt, userinfo;
  int window,  offset,   graph,  resx,    resy;
  int gfrom,   gfromSave,gto,    gtoSave, golap, glen;
  char printopt;
  char outputname[FILENAME_MAX+1], parname[FILENAME_MAX+1];
  char blastArg[FILENAME_MAX+1],   grnameArg[FILENAME_MAX+1];
  char matname[FILENAME_MAX+1];
  char grname[FILENAME_MAX+1],     fstname[FILENAME_MAX+1];

  // Par File
  FILE *fp;
  char   clef[20];
  char   versionPAR[FILENAME_MAX+1];
  double FsP,StartP,StartB,StopP,TransStartP,TransStopP;
  double AccP[2],AccB[2],DonP[2],DonB[2],BlastS[8],EstP;
  double TransitionWeight[5];
  int EstM;
  int MinFivePrime, MinThreePrime;
  int MinEx, MinIn, MinSg, MinFlow, MinConv, MinDiv;
  int MinLength[18];  // Les longueurs

  Parameters  ();
  ~Parameters ();
  void InitParam (int, char *[]);
};

#endif
