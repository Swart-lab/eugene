//=======================================================================
//            Copyright (c) 2003 by INRA. All rights reserved            
//                Redistribution is not permitted without                
//                the express written permission of INRA                 
//                   eMail : tschiex@toulouse.inra.fr                    
//------------------------------------------------------------------------
// File        : EuGeneTk/Parametrization/LineSearch/Chromosome.cc             
// Description : Class Data and Chromosome of the Genetic algorithm
// Authors     : P.Bardou, S.Foissac, M.J.Cros, A.Moisan, T.Schiex       
//=======================================================================



#include "Chromosome.h"

#include "../ParaOptimization.h"
#include "Genetic.h"

extern ParaOptimization OPTIM;


//-------------------------------------------------------
//-------------------------------------------------------
Data::Data(int n)
{
  Genetic* algo = (Genetic*) OPTIM.Algorithms[OPTIM.AlgoIndex];

  for (int i=0; i<n; i++)
    P.push_back( algo->ParMin[i] + 
		 algo->Rand->RandUniform()*(algo->ParMax[i] - algo->ParMin[i]) );
}




//-------------------------------------------------------
// Data distance for clustering
//-------------------------------------------------------
double Data::GetDistance(Data* d)
{
  Genetic* algo = (Genetic*) OPTIM.Algorithms[OPTIM.AlgoIndex];

  double dist=0.0;

  for (unsigned int i=0; i< algo->ParName.size(); i++)
    dist += ( (algo->ParMax[i]==algo->ParMin[i])
	      ? 0 
	      : pow( (d->P[i] - P[i])/(algo->ParMax[i] - algo->ParMin[i]), 2)  );
 
  return dist;
}



//-------------------------------------------------------
//-------------------------------------------------------
void Data::Recenter(void)
{
  Genetic* algo = (Genetic*) OPTIM.Algorithms[OPTIM.AlgoIndex];

  for (unsigned int i=0; i<P.size(); i++) {
    if (P[i] < algo->ParMin[i]) P[i] = algo->ParMin[i];
    if (P[i] > algo->ParMax[i]) P[i] = algo->ParMax[i];
  }
}


//-------------------------------------------------------
//-------------------------------------------------------
void Data::PutBarycenter(int nb1, Data *data2, int nb2) 
{
  if((nb1+nb2)==0){
    nb1=1;nb2=1;
  }

  for (unsigned int i=0; i<P.size(); i++) 
    P[i] = (P[i]*nb1 + data2->P[i]*nb2)/(nb1+nb2);

  Recenter();
}


//-------------------------------------------------------
// Constructor
//-------------------------------------------------------
Chromosome::Chromosome(int n)
{
  RawFitness = 0;
  ScaledFitness = 0;
  IsEvaluated = false;
  IsNew = true;
  IsBestInCluster = false;
  NoCluster = 0;
  data = new Data(n);
}

//-------------------------------------------------------
// Destructor
//-------------------------------------------------------
Chromosome::~Chromosome(void)
{
  delete data;
}


//-------------------------------------------------------
// Evaluate chromosome
//-------------------------------------------------------
void Chromosome::Evaluate(void)
{
  Genetic* algo = (Genetic*) OPTIM.Algorithms[OPTIM.AlgoIndex];
  
  algo->Par = data->P;
  for (unsigned int k=0; k<algo->Para.size(); k++) algo->Para[k] = algo->Par[algo->ParaPar[k]];
  RawFitness = OPTIM.ParaEvaluate();

  IsEvaluated = true;
}


//-------------------------------------------------------
// Constructor
//-------------------------------------------------------
Cluster::Cluster(int n)
{
  CentralData = new Data(n);
}


//-------------------------------------------------------
// Destructor
//-------------------------------------------------------
Cluster::~Cluster(void)
{
  delete CentralData;
}
