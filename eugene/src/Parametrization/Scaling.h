//=======================================================================
//            Copyright (c) 2003 by INRA. All rights reserved            
//                Redistribution is not permitted without                
//                the express written permission of INRA                 
//                   eMail : tschiex@toulouse.inra.fr                    
//------------------------------------------------------------------------
// File        : EuGeneTk/Parametrization/LineSearch/Scaling.h             
// Description : Class realizing the scaling for the genetic algorithm
// Authors     : P.Bardou, S.Foissac, M.J.Cros, A.Moisan, T.Schiex       
//=======================================================================

#ifndef SCALING_H_INCLUDED
#define SCALING_H_INCLUDED

#include "../OptiAlgorithm.h"


class Scaling {

 public:
  int Type;

  Scaling(int type);

  void Scalestat(int gen);
};

#endif
