//=======================================================================
//            Copyright (c) 2003 by INRA. All rights reserved            
//                Redistribution is not permitted without                
//                the express written permission of INRA                 
//                   eMail : tschiex@toulouse.inra.fr                    
//------------------------------------------------------------------------
// File        : EuGeneTk/Parametrization/OptiAlgorithm.h             
// Description : the mother class of optimization algorithms classses
// Authors     : P.Bardou, S.Foissac, M.J.Cros, A.Moisan, T.Schiex       
//=======================================================================

#ifndef OPTI_ALGORITHM_H_INCLUDED
#define OPTI_ALGORITHM_H_INCLUDED


#include <vector>
#include <string>

#include "Random.h"


// Definition of cluster relations
enum RELATION {LINKED, IDENTICAL, NO_RELATION};
typedef enum RELATION CLUSTER_RELATION;


class OptiAlgorithm {

 public:
  vector <double> Para;  
  vector <string> ParaName;
  vector <double> ParaMax;
  vector <double> ParaMin; 
  vector <int> ParaCluster;
  vector < vector <int> > ParaClusters;
  vector <CLUSTER_RELATION> ParaClusterRelations;
  int NbParaCluster;
  bool IsTracing;
  Random* Rand;

  OptiAlgorithm(void);
  virtual ~OptiAlgorithm(void) {};
  virtual void Optimize(bool is_chaining) = 0;
  string ReduceName(string s);
};




#endif
