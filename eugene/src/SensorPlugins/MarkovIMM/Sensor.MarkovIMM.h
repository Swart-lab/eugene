#ifndef  SENSOR_MARKOVIMM_H_INCLUDED
#define  SENSOR_MARKOVIMM_H_INCLUDED

#include "../../EuGene/Sensor.h"
#include "../0_SensorTk/BStrArray.h"

/****************************************************************
 **                     SensorMarkovIMM                        **
 ****************************************************************/
class SensorMarkovIMM : public Sensor
{
 private:
  // to avoid to reread the IMM matrices, they are stored in the list IMMatrixList
  // associated to a matrix: 
  //    - the interval of GC% [minGC maxGC] where it is defined
  //    - the request of using model M0 for IG
  // for a given instance, IMMatrix_index allows to know in lists which matrice 
  // and relative information to consider
  static std::vector < std::vector <BString_Array*> > IMMatrixList;  
  static std::vector <double> minGCList;
  static std::vector <double> maxGCList;
  static std::vector <bool>   UseM0asIGList;
  const static int  MODEL_LEN = 9;
  const static int  SIMPLE_MODEL_LEN = 6;
  const static int  ALPHABET_SIZE = 4;
  int IMMatrix_index;
  

 public:
  SensorMarkovIMM         (int n, DNASeq *X);
  virtual ~SensorMarkovIMM();
  virtual void Init       (DNASeq *);
  virtual void GiveInfo   (DNASeq *, int, DATA *);
  virtual void Plot       (DNASeq *);
  virtual void PostAnalyse(Prediction *);
};

extern "C" SensorMarkovIMM * builder0( int n, DNASeq *X) { return new SensorMarkovIMM(n, X);}

// reserve memory for static variables
std::vector < std::vector <BString_Array*> > SensorMarkovIMM::IMMatrixList;
std::vector <double>                         SensorMarkovIMM::minGCList;
std::vector <double>                         SensorMarkovIMM::maxGCList;
std::vector <bool>                           SensorMarkovIMM::UseM0asIGList;

#endif
