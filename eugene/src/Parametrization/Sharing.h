//=======================================================================
//            Copyright (c) 2003 by INRA. All rights reserved            
//                Redistribution is not permitted without                
//                the express written permission of INRA                 
//                   eMail : tschiex@toulouse.inra.fr                    
//------------------------------------------------------------------------
// File        : EuGeneTk/Parametrization/LineSearch/Sharing.h             
// Description : Class realizing the sharing for the genetic algorithm
// Authors     : P.Bardou, S.Foissac, M.J.Cros, A.Moisan, T.Schiex       
//=======================================================================

#ifndef SHARING_H_INCLUDED
#define SHARING_H_INCLUDED


#include "../OptiAlgorithm.h"


class Sharing {

 public:
  double Value;

  Sharing(double value);
  void Share (void);

 private:
  double EvalShare(double distance);

};

#endif
