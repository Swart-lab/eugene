//=======================================================================
//            Copyright (c) 2003 by INRA. All rights reserved            
//                Redistribution is not permitted without                
//                the express written permission of INRA                 
//                   eMail : tschiex@toulouse.inra.fr                    
//------------------------------------------------------------------------
// File        : EuGeneTk/Parametrization/LineSearch/Mutating.cc         
// Description : Class realizing mutation for the genetic algorithm
// Authors     : P.Bardou, S.Foissac, M.J.Cros, A.Moisan, T.Schiex       
//=======================================================================



#include "Mutating.h"

#include "../ParaOptimization.h"
#include "Genetic.h"

extern ParaOptimization OPTIM;


//-------------------------------------------------------
// Constructor
//-------------------------------------------------------
Mutating::Mutating(double proba)
{
  Chromosome* c;
  Genetic* algo =  (Genetic*) OPTIM.Algorithms[OPTIM.AlgoIndex];

  Proba = proba;
  if (algo->IsSAMutating) 
      for (int i=0; i<2; i++) {
	c = new Chromosome(algo->ParName.size());
	ChildrenPopulation.push_back(c);
      }
}


//-------------------------------------------------------
// Destructor
//-------------------------------------------------------
Mutating::~Mutating(void)
{
  Genetic* algo =  (Genetic*) OPTIM.Algorithms[OPTIM.AlgoIndex];

  if (algo->IsSAMutating) 
    for (int i=0; i<2; i++) 
      delete ChildrenPopulation[i];
}





//-------------------------------------------------------
//-------------------------------------------------------
void Mutating::Muteval (void)
{
  Genetic* algo =  (Genetic*) OPTIM.Algorithms[OPTIM.AlgoIndex];

  int i,j,n,k;

  n = (int)(algo->NbElement*Proba + 0.5);  

  if (!algo->IsSAMutating)
    for (i=0; i<n; i++) {
      ChooseChromosome(&j);
      Mutation(algo->Population[j]) ;
      algo->Population[j]->IsNew = true;
      algo->Population[j]->IsEvaluated = false;
    }
  else {
    for (i=0; i<n; i++){
      ChooseChromosome(&j) ;
      ChildrenPopulation[0]->data->P = algo->Population[j]->data->P;
      Mutation(ChildrenPopulation[0]) ;
      algo->Par = ChildrenPopulation[0]->data->P;
      for (k=0; k<(int)algo->Para.size(); k++) algo->Para[k] = algo->Par[algo->ParaPar[k]];
      ChildrenPopulation[0]->RawFitness = OPTIM.ParaEvaluate();
      /* Here everything as been calculated, we just have to crossover */
      algo->SA->SATournament2(ChildrenPopulation[0], algo->Population[j]);
      algo->Population[j]->IsEvaluated = true;
    }
  }
} 



//-------------------------------------------------------
//-------------------------------------------------------
void Mutating::ChooseChromosome(int* n)
{
  Genetic* algo =  (Genetic*) OPTIM.Algorithms[OPTIM.AlgoIndex];

  int k;

  /* If we use an elitist strategy, we can not mutate the overall */
  /* best chromosom (we keep at least one best overall) */
  /* Anyway, we do not mutate an element which has just been added */
  /* in the population */
  if (algo->IsElitist) {
    k=0;
    do
      *n = (int)(algo->Rand->RandUniform()*(double)algo->NbElement);
    while ( ((algo->Population[*n] == algo->BestChromosome)
	     ||(algo->Population[*n]->IsNew)
	     ||(algo->Population[*n]->IsBestInCluster) )
	    && (k++<1000) );
    /* If we went over limit then we must choose one but we */
    /* refuse best_chrom, we accept best_of_cluster to be   */
    /* choosen and then deleted */
    if (k>1000) {
      do
	*n = (int)(algo->Rand->RandUniform()*(double)algo->NbElement);
      while ( (algo->Population[*n] == algo->BestChromosome)
	      ||(algo->Population[*n]->IsNew) ) ;
    }
  } else
    do
      *n = (int)(algo->Rand->RandUniform()*algo->NbElement);
    while (algo->Population[*n]->IsNew);
}


//-------------------------------------------------------
//-------------------------------------------------------
void Mutating::Mutation(Chromosome* c)
{
  Genetic* algo =  (Genetic*) OPTIM.Algorithms[OPTIM.AlgoIndex];

  double alpha=0.2;
  // 68% des valeurs mutees tombent dans l'intervalle 
  // [ Valeur initiale-alpha*(Max-Min) -> Valeur initiale+alpha*(Max-Min) ]
  // et 95% entre Val.init-2alpha*(Max-Min) et Val.init+2alpha*(Max-Min)
  // amelioration: diminuer alpha au cours du temps?
  
  /* Evalue la donnée bruitée avec un bruit gaussien*/
  for (unsigned int i=0; i<c->data->P.size(); i++)
    c->data->P[i] = c->data->P[i] + 
      algo->Rand->RandGauss()*alpha*(algo->ParMax[i] - algo->ParMin[i]);

  c->data->Recenter();
}
