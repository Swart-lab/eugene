//=======================================================================
//            Copyright (c) 2003 by INRA. All rights reserved            
//                Redistribution is not permitted without                
//                the express written permission of INRA                 
//                   eMail : tschiex@toulouse.inra.fr                    
//------------------------------------------------------------------------
// File        : EuGeneTk/Parametrization/LineSearch/SimulatedAnnealing.h             
// Description : Class realizing the simulated annealing for the genetic algorithm
// Authors     : P.Bardou, S.Foissac, M.J.Cros, A.Moisan, T.Schiex       
//=======================================================================

#ifndef SIMULTED_ANNEALING_H_INCLUDED
#define SIMULTED_ANNEALING_H_INCLUDED


#include "../Random.h"
#include "Chromosome.h"


class SimulatedAnnealing {

 public:
  SimulatedAnnealing(void);
  void SAInit(void);
  void SAUpdate(int gen);
  void SATournament4(Chromosome *C1,Chromosome *C2, 
		     Chromosome *P1, Chromosome *P2);
  void SATournament2(Chromosome* c1, Chromosome* c2);

 private:
  int CP;
  double CC1;
  double CC2;
  double CC;  /* CC : weight value for temperature modification during */
              /*      the annealing scheme : T(k+1)=CC*T(k) where k    */
              /*      is the current generation number. Two annealing  */
              /*      schemes are used during the temperature decrease:*/
              /*      one with CC1 and the other with CC2.             */
  double Tx;  /* we change the annealing scheme when we reach */
              /* the temperature Tx                           */
  double T; /* re-annealing temperature */
};

#endif
