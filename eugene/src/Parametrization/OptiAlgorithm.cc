//=======================================================================
//            Copyright (c) 2003 by INRA. All rights reserved            
//                Redistribution is not permitted without                
//                the express written permission of INRA                 
//                   eMail : tschiex@toulouse.inra.fr                    
//------------------------------------------------------------------------
// File        : EuGeneTk/Parametrization/OptiAlgorithm.cc       
// Description : the mother class of optimization algorithms classses
// Authors     : P.Bardou, S.Foissac, M.J.Cros, A.Moisan, T.Schiex       
//=======================================================================

#include <iostream>
#include <iomanip>


#include "OptiAlgorithm.h"


#include "../EuGene/Param.h"

extern Parameters PAR;

//-------------------------------------------------------
// Constructor
//-------------------------------------------------------
OptiAlgorithm::OptiAlgorithm(void)
{
  int n, c;
  std::string relation;
  std::vector <int> v;

  IsTracing = ( (((std::string) PAR.getC("ParaOptimization.Trace")) == "TRUE") ? true : false );

  NbParaCluster = PAR.getI("ParaOptimization.NbCluster");
  for (int i = 0; i < NbParaCluster; i++) {
    ParaClusters.push_back( v );
    relation = (std::string) PAR.getC("ParaOptimization.Cluster", i);
    if (relation == "IDENTICAL")
      ParaClusterRelations.push_back(IDENTICAL); 
    else
      if (relation == "LINKED")
	ParaClusterRelations.push_back(LINKED); 
      else 	{
	std::cerr <<"ERROR: Bad cluster relation "<<relation<<"  in the parameter file."<<std::endl; 
	exit(100);
      }
  }

  n = PAR.getI("ParaOptimization.NbParameter");
  for (int i=0; i<n; i++) {
    ParaName.push_back( PAR.getC("ParaOptimization.Para.Name", i) );
    ParaMin.push_back( PAR.getD("ParaOptimization.Para.Min" , i) );
    ParaMax.push_back( PAR.getD("ParaOptimization.Para.Max" , i) );
    c = PAR.getI("ParaOptimization.Para.Cluster", i);
    ParaCluster.push_back( c );
    ParaClusters[c].push_back( i );
  }
  for (int i=0; i<n; i++) 
    if ( ParaName[i][ParaName[i].size()-1] != '*' ) {
      std::cerr <<"ERROR: Parameters to optimize must have a name finishing by '*'."<<std::endl;
      exit(100);
    }

  // Set numerical precision for displayed values
  std::cout.setf(std::ios::showpoint);
  std::cout << setiosflags (std::ios::fixed);
  std::cout.precision(4);
}

//-------------------------------------------------------
// return a string with maximum 6 signs
// if more than 6 letters returns: the 2 first letters + '..' + the 2 last
//-------------------------------------------------------
std::string OptiAlgorithm::ReduceName(std::string s)
{
  std::string ss, s0, s1, sn2, sn1;
  int n = s.length();
  
  if (n > 6) {
   s0 = s[0]; s1 = s[1]; sn2 = s[n-2]; sn1 = s[n-1];
   ss = s0 + s1 + ".."+ sn2 + sn1;
 } else
   ss =s;
 
 return ss;
}

