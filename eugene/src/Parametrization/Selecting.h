//=======================================================================
//            Copyright (c) 2003 by INRA. All rights reserved            
//                Redistribution is not permitted without                
//                the express written permission of INRA                 
//                   eMail : tschiex@toulouse.inra.fr                    
//------------------------------------------------------------------------
// File        : EuGeneTk/Parametrization/LineSearch/Selecting.h             
// Description : Class realizing the selection for the genetic algorithm
// Authors     : P.Bardou, S.Foissac, M.J.Cros, A.Moisan, T.Schiex       
//=======================================================================

#ifndef SELECTING_H_INCLUDED
#define SELECTING_H_INCLUDED




class Selecting {

 public:
  int Type;

  Selecting(int type);
  ~Selecting(void);
  void Selection(void);

 private:
  double* tab;
  double* fraction;
  int* choices;
  int nremain;

};

#endif
