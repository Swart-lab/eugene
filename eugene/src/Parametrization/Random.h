//=======================================================================
//            Copyright (c) 2003 by INRA. All rights reserved            
//                Redistribution is not permitted without                
//                the express written permission of INRA                 
//                   eMail : tschiex@toulouse.inra.fr                    
//------------------------------------------------------------------------
// File        : EuGeneTk/Parametrization/Random.h             
// Description : Class for generators of numbers
// Authors     : P.Bardou, S.Foissac, M.J.Cros, A.Moisan, T.Schiex       
//=======================================================================

#ifndef RANDOM_H_INCLUDED
#define RANDOM_H_INCLUDED


#include <string>


class Random {

 public:
  string Seed;

  Random(string seed);
  Random(void);
  int RandBoundedInt (int bound);
  double RandUniform(void);
  double RandGauss(void);
  int Flip(double prob);

 private:
  bool IsInitialized;
  unsigned short int xsubi[3];
  int iset;
  double gset;

};

#endif
