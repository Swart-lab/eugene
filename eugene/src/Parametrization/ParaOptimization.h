//=======================================================================
//            Copyright (c) 2003 by INRA. All rights reserved            
//                Redistribution is not permitted without                
//                the express written permission of INRA                 
//                   eMail : tschiex@toulouse.inra.fr                    
//------------------------------------------------------------------------
// File        : EuGeneTk/Parametrization/ParaOptimization.h             
// Description : The ParaOptimization class optimizes EuGene parameters
// Authors     : P.Bardou, S.Foissac, M.J.Cros, A.Moisan, T.Schiex       
//=======================================================================

#ifndef PARA_OPTIMIZATION_H_INCLUDED
#define PARA_OPTIMIZATION_H_INCLUDED


#include <vector>
#include <string>

#include "OptiAlgorithm.h"
#include "../EuGene/MSensor.h"
#include "../EuGene/DNASeq.h"

class ParaOptimization {

 public:
  std::vector <OptiAlgorithm*> Algorithms;
  int                          AlgoIndex;     // index of running algorithm

  ~ParaOptimization(void);
  void ParaOptimize (int argc, char* argv[]);
  double ParaEvaluate (void);  

 private:
  std::vector <MasterSensor*>  MSensors;
  std::vector <DNASeq*>        Sequences;
  std::vector<std::string>     SeqNames;
  std::string                  TrueCoordFile;
  std::string                  ExecutableName; // first part of parameter file
  bool                         IsTest;         // if true then test mode: 
                                      // do not read sequences and related info
                                      // evaluate para with the method NormalLaw

  void Init(int argc, char* argv[]);
  double NormalLaw(double x);
};


#endif
