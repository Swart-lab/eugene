//=======================================================================
//            Copyright (c) 2003 by INRA. All rights reserved            
//                Redistribution is not permitted without                
//                the express written permission of INRA                 
//                   eMail : tschiex@toulouse.inra.fr                    
//------------------------------------------------------------------------
// File        : EuGeneTk/Parametrization/LineSearch/CrossingOver.h             
// Description : Class for cross over of the Genetic algorithm
// Authors     : P.Bardou, S.Foissac, M.J.Cros, A.Moisan, T.Schiex       
//=======================================================================

#ifndef CROSSOVER_H_INCLUDED
#define CROSSOVER_H_INCLUDED

#include <vector>


#include "Chromosome.h"


class CrossingOver {
 public:
  double Proba;

  CrossingOver(double proba);
  ~CrossingOver(void);
  void Crosseval(void); 

 private:
  vector<Chromosome*> ChildrenPopulation;

  void Crossover(Chromosome* c1, Chromosome* c2);
  void ChooseCouple(int* n1, int* n2);

};

#endif
