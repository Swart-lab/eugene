//=======================================================================
//            Copyright (c) 2003 by INRA. All rights reserved            
//                Redistribution is not permitted without                
//                the express written permission of INRA                 
//                   eMail : tschiex@toulouse.inra.fr                    
//------------------------------------------------------------------------
// File        : EuGeneTk/Parametrization/LineSearch/Mutating.h             
// Description : Class realizing mutation for the genetic algorithm
// Authors     : P.Bardou, S.Foissac, M.J.Cros, A.Moisan, T.Schiex       
//=======================================================================

#ifndef MUTATING_H_INCLUDED
#define MUTATING_H_INCLUDED


#include <vector>
#include <values.h>

#include "Chromosome.h"


class Mutating {

 public:
  double Proba;

  Mutating(double proba);
  ~Mutating(void);
  void Muteval (void); 

 private:
  vector<Chromosome*> ChildrenPopulation;
  void ChooseChromosome(int* n);
  void Mutation(Chromosome* c);
};

#endif
