#ifndef WAM_INCLUDED
#define WAM_INCLUDED
//#include "../../EuGene/Sensor.h"
#include "./markov.h"

#ifndef FILENAME_MAX
#define FILENAME_MAX        1024
#endif

#define TPFILESUFFIX ".TP."
#define FPFILESUFFIX ".FP."
#define SUFFIXLENGTH 4

/*************************************************************
 **                 WAM                    **
 *************************************************************/

class WAM
{
 private:
  int MarkovianOrder;
  int MotifLength;
  Chaine* Alphabet;
  TabChaine<Chaine, unsigned short int>* TPMOD; // two tab. of markovian models
  TabChaine<Chaine, unsigned short int>* FPMOD;

 public:
  WAM();
  WAM(int order, int length, char* alphabet, char* prefixfilename);
  ~WAM();
  double ScoreTheMotif(char* oligont);  // sum of likelihood ratio for each position of the entire motif
};

#endif
