//=======================================================================
//            Copyright (c) 2003 by INRA. All rights reserved            
//                Redistribution is not permitted without                
//                the express written permission of INRA                 
//                   eMail : tschiex@toulouse.inra.fr                    
//------------------------------------------------------------------------
// File        : EuGeneTk/Parametrization/LineSearch.cc           
// Description : Class of the Line Search algorithm
// Authors     : P.Bardou, S.Foissac, M.J.Cros, A.Moisan, T.Schiex       
//=======================================================================

#include <string>

#include "LineSearch.h"


#include "../../EuGene/Param.h"
#include "../ParaOptimization.h"
#include "../Random.h"

extern Parameters PAR;
extern ParaOptimization OPTIM;

//-------------------------------------------------------
// Constructor
//-------------------------------------------------------
LineSearch::LineSearch (void) : OptiAlgorithm()
{
  double step;
  int i;
  vector <bool> FirstEltCluster;

  NbMaxCycle = PAR.getI("LineSearch.NbMaxCycle");
  NbMinCycle = PAR.getI("LineSearch.NbMinCycle");
  NbMaxStab = PAR.getI("LineSearch.NbMaxStab");
  DivInter = PAR.getI("LineSearch.DivInter");
  Alpha = PAR.getD("LineSearch.Alpha");
  EvolutionMini = PAR.getD("LineSearch.EvolutionMini");
  
  for (int i=0; i<PAR.getI("ParaOptimization.NbParameter"); i++) {
    Para.push_back( PAR.getD("LineSearch.Para.Init", i) );
    ParaMinStep.push_back( PAR.getD("LineSearch.Para.Step", i) );
    ParaMinInter.push_back( PAR.getD("LineSearch.Para.MinInit", i) );
    ParaMaxInter.push_back( PAR.getD("LineSearch.Para.MaxInit", i) );
    ParaLInter.push_back( ParaMaxInter[i] - ParaMinInter[i] );
  }
    
  for (unsigned int i=0; i< Para.size(); i++) {
    step = (ParaMaxInter[i] - ParaMinInter[i]) / DivInter;
    ParaStep.push_back( (step > ParaMinStep[i]) ? step : ParaMinStep[i]);
  }
  
  Rand = new Random(PAR.getC("LineSearch.Seed"));

  // Update de names of parameters with just the first elt of clusters IDENTICAL
  for (i=0; i<(int)ParaClusters.size(); i++) FirstEltCluster.push_back( true );
  MsgParaNames="";
  for (i=0; i<(int)ParaName.size(); i++) 
    if (ParaClusterRelations[ParaCluster[i]] == IDENTICAL) {
      if ( FirstEltCluster[ParaCluster[i]] ) {
	MsgParaNames += ReduceName(ParaName[i]) + "\t";
	FirstEltCluster[ParaCluster[i]] = false;
      }
    } else
      MsgParaNames += ReduceName(ParaName[i]) + "\t";

}



/*======================================================================*/
/* Optimize: launches optimization                                      */
/*======================================================================*/
void LineSearch :: Optimize(bool is_chaining) 
{
  int n;
  int STOP;
  double FitOptPrec = Fitness;
  FitOpt = 0.0;
  unsigned int i; int c;
  string warning_message;

  // if chains another algo then takes the optimum found as initial point
  if (is_chaining)
    Para = OPTIM.Algorithms[OPTIM.AlgoIndex-1]->Para;

  // Display parametrization of the algorithm
  cout <<endl << "---------------------------------------------------------------"<<endl;
  cout << "Optimization of EuGène parameters with the LineSearch algorithm"<<endl<<endl;
  cout << "---------------------------------------------------------------"<<endl;
  cout << "Parametrisation of the algorithm:"<<endl<<endl;
  cout << "NbMaxCycle: " << NbMaxCycle << "\tNbMinCycle: " << NbMinCycle 
       << "\tNbMaxStab: " << NbMaxStab <<endl;
  cout << "DivInter: " << DivInter << "\tAlpha: " << Alpha 
       << "\tEvolutionMini: " << EvolutionMini <<endl;
  cout << "Seed: " << Rand->Seed <<endl;
  cout << "Trace: " << ((IsTracing) ? "TRUE" : "FALSE") <<endl<<endl;

  cout << "NbParameter: " << (int) Para.size() << "\tNbCluster: " << NbParaCluster <<endl;
  cout << "Param: \t"; for (i=0; i<ParaName.size(); i++) cout << ReduceName(ParaName[i]) << "\t"; cout <<endl;
  cout << "Init: \t"; for (i=0; i<Para.size(); i++) cout << Para[i] << "\t"; cout <<endl;
  cout << "Min: \t"; for (i=0; i<ParaMin.size(); i++) cout << ParaMin[i] << "\t"; cout <<endl;
  cout << "Max: \t"; for (i=0; i<ParaMax.size(); i++) cout << ParaMax[i] << "\t"; cout <<endl;
  cout << "MinStep: \t"; for (i=0; i<ParaMinStep.size(); i++) cout << ParaMinStep[i] << "\t"; cout <<endl;

  cout <<endl;
  for (c=0; c<NbParaCluster; c++) {
    cout << "Cluster[" << c << "]: \t";
    for (unsigned int j=0; j<ParaClusters[c].size(); j++) 
      cout << ParaName[((ParaClusters[c])[j])] <<" ";
    if (ParaClusterRelations[c] == IDENTICAL) 
      cout << "   IDENTICAL"<<endl;
    else
      cout << "   LINKED"<<endl;
  }

  Fitness = OPTIM.ParaEvaluate();
  cout <<endl<< "Fitness of initial point: " << Fitness <<endl<<endl;

  // general loop
  n=0; STOP = 0;
  while ((n < NbMaxCycle) && ((n < NbMinCycle) || (STOP < NbMaxStab))) {
    cout <<endl <<"Cycle: " << (n + 1)
	 << " -----------------------------------------------------\n";
    cout << MsgParaNames << endl;
    
    for (int k = 0; k < NbParaCluster; k++) { 
      ScanCluster(k);
      ChooseOptimal();

      if (IsTracing) { cout << "Local optimal point:" <<endl; PrintParam(); cout <<endl;}
    }
    
    if ((FitOpt - FitOptPrec) < EvolutionMini)
      STOP += 1; 
    else { 
      STOP = 0;
      FitOptPrec = FitOpt; 
    }
     
    ReduceSearch();
    n++;
  }
  
  cout <<endl<< "---------------------------------------------------------------"<<endl;
  cout << "LineSearch stops after " << n << " cycles."<<endl;
  if (n==NbMaxCycle) cout <<"Maximum number of cycles achieved."<<endl;
  if (STOP==NbMaxStab) cout <<"Fitness is stable since "<<NbMaxStab<<" cycles." <<endl;
  cout <<endl << "Final Optimal Point:" <<endl; 
  for (i=0; i<ParaName.size(); i++) cout << ReduceName(ParaName[i]) << "\t"; cout <<endl;
  for (i=0; i<ParaName.size(); i++) cout << Para[i] << "\t"; cout <<endl;
  cout << "---------------------------------------------------------------"<<endl;
  cerr <<endl<< warning_message ;
}


/*======================================================================*/
/* UpdateOpt: Update list of points which fitness is equal to current   */
/*            optimal fitness                                           */
/*======================================================================*/ 
void LineSearch :: UpdateOpt(void) 
{  
  if (FitOpt <= Fitness) { 
    if (FitOpt < Fitness) {
      FitOpt = Fitness;
      Optimums.clear(); 
    }
    
    Optimums.push_back(Para); 
  }
}


/*======================================================================*/
/* ScanCluster: scans parameters of current cluster                     */
/*======================================================================*/
void LineSearch :: ScanCluster(int k) 
{  
  int n, m;
  unsigned int i;

  if (ParaClusterRelations[k] == IDENTICAL) { 
    // We assume that identical parameters do have identical intervals
    n = (ParaClusters[k])[0];
	  
    for (double t = ParaMinInter[n]; t <= ParaMaxInter[n]; t += ParaStep[n]) {
      // all the parameters of the cluster are set the same value
      for (i=0; i<ParaClusters[k].size(); i++) {
	m = (ParaClusters[k])[i];
	Para[m] = t;
      }
      
      Fitness = OPTIM.ParaEvaluate();
      if (IsTracing) PrintParam();
      UpdateOpt(); 
    }
  } else
    CartesianProduct(0,k);
}


/*======================================================================*/
/* CartesianProduct: used in ScanCluster(k) when cluster[k] is linked   */
/*======================================================================*/
void LineSearch :: CartesianProduct(int j, int k) 
{
  int m = (ParaClusters[k])[j];
  
  if ( j == ((int) ParaClusters[k].size()) - 1 ) 
    for (double s = ParaMinInter[m]; s <= ParaMaxInter[m]; s += ParaStep[m]) {
      Para[m] = s;
      Fitness = OPTIM.ParaEvaluate();
      UpdateOpt();
      if (IsTracing) PrintParam();
    }
  else {
    j++;
    for (double t = ParaMinInter[m]; t <= ParaMaxInter[m]; t += ParaStep[m]) {
      Para[m] = t;
      CartesianProduct(j,k);
    }
  }
}

	  
/*======================================================================*/
/* ChooseOptimal: chooses randomly one optimal point from Optimums vector    */
/*======================================================================*/
void LineSearch :: ChooseOptimal(void) 
{  
  int n = (int) Rand->RandBoundedInt((int) Optimums.size());
  
  for (unsigned int i=0; i<Para.size(); i++) Para[i] = (Optimums[n])[i];
  
  Optimums.clear();
  Optimums.push_back(Para);
  Fitness = FitOpt;
}


/*======================================================================*/
/* ReduceSearch: Reduces interval search for each parameter             */
/*======================================================================*/
void LineSearch :: ReduceSearch(void) {
  
  double  marge, step; 
  
  for (unsigned int i=0; i<Para.size(); i++) {
    // Reduce interval if possible
    ParaLInter[i] *= Alpha;
    marge = ParaLInter[i]/2;    
    if (ParaLInter[i] > ParaMinStep[i]) {
      ParaMinInter[i] = ((ParaMin[i]>(Para[i] - marge)) ? ParaMin[i] : (Para[i]-marge));
      ParaMaxInter[i] = ((ParaMax[i]<(Para[i] + marge)) ? ParaMax[i] : (Para[i]+marge));
    }
    // Update the step to explore the interval

    step = (ParaMaxInter[i] - ParaMinInter[i]) / DivInter;
    ParaStep[i] = ( (step > ParaMinStep[i]) ? step : ParaMinStep[i] );
  }
}



//-------------------------------------------------------
// Print Para values and associated fitness
// BEWARE: just the first para of an IDENTICAL cluster is printed.
//-------------------------------------------------------
void LineSearch::PrintParam(void)
{
  unsigned int i;
  vector <bool> FirstEltCluster;

  // parameters value with just the first elt of clusters IDENTICAL
  for (i=0; i<ParaClusters.size(); i++) FirstEltCluster.push_back( true );
  for (i=0; i<ParaName.size(); i++) 
    if (ParaClusterRelations[ParaCluster[i]] == IDENTICAL) {
      if ( FirstEltCluster[ParaCluster[i]] ) {
	cout << Para[i] << "\t";
	FirstEltCluster[ParaCluster[i]] = false;
      }
    } else
      cout << Para[i] << "\t";

  cout << "Fitness = " << Fitness << "\n";
}

