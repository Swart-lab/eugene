//=======================================================================
//            Copyright (c) 2003 by INRA. All rights reserved            
//                Redistribution is not permitted without                
//                the express written permission of INRA                 
//                   eMail : tschiex@toulouse.inra.fr                    
//------------------------------------------------------------------------
// File        : EuGeneTk/Parametrization/LineSearch/LineSearch.h             
// Description : Class of the Line Search algorithm
// Authors     : P.Bardou, S.Foissac, M.J.Cros, A.Moisan, T.Schiex       
//=======================================================================

#ifndef LINESEARCH_H_INCLUDED
#define LINESEARCH_H_INCLUDED


#include "../OptiAlgorithm.h"


class LineSearch : public OptiAlgorithm {

 private:
  std::vector <double> ParaMinInter;
  std::vector <double> ParaMaxInter;
  std::vector <double> ParaLInter;
  std::vector <double> ParaStep;
  std::vector <double> ParaMinStep;
  std::vector < std::vector <double> > Optimums;
  int NbMaxCycle;
  int NbMinCycle;
  int NbMaxStab;
  int DivInter;
  double Alpha;
  double EvolutionMini;
  double FitOpt;
  double Fitness;
  std::string MsgParaNames;

  void ReduceSearch(void);
  void ScanCluster(int k);
  void CartesianProduct(int j, int k);
  void ChooseOptimal(void);
  void UpdateOpt(void);
  void PrintParam(void);
  
 public:

  LineSearch(void);

  void Optimize(bool is_chaining);
};


#endif
