//=======================================================================
//            Copyright (c) 2003 by INRA. All rights reserved            
//                Redistribution is not permitted without                
//                the express written permission of INRA                 
//                   eMail : tschiex@toulouse.inra.fr                    
//------------------------------------------------------------------------
// File        : EuGeneTk/Parametrization/LineSearch/Genetic.h             
// Description : Class of the Genetic algorithm
// Authors     : P.Bardou, S.Foissac, M.J.Cros, A.Moisan, T.Schiex       
//=======================================================================

#ifndef GENETIC_H_INCLUDED
#define GENETIC_H_INCLUDED


#include "../OptiAlgorithm.h"

#include "CrossingOver.h"
#include "Mutating.h"
#include "Selecting.h"
#include "Scaling.h"
#include "Sharing.h"
#include "SimulatedAnnealing.h"
#include "Chromosome.h"


class Genetic : public OptiAlgorithm {

 private:
  bool IsScaling;

  int NbRun;
  double FitnessOpti;
  std::vector<double> ParOpti;

  CrossingOver* CrossOver;
  Mutating* Mutate;
  Selecting* Select;
  Scaling* Scale;
  Sharing* Share;

  std::vector<int> IndexesOfBestInClusters;
  double d_mean;
  double minmaxfactor;

  void InitPopulation(void);
  void EvalPopulation(void);
  int MakeClusters(void);
  void MarkBestOfClusters(void);
  int MergeClusters(int cl1, int cl2,int last);
  int CountNopt(int nb_clusters, double max);

 public:
  int NbElement;
  int NbGeneration;
  int NbCluster;       // Nb of effective AG clusters
  bool IsElitist;
  bool IsClustering;
  bool IsSharing;
  bool IsSACrossingOver;
  bool IsSAMutating;
  double Elitism;
  Chromosome* BestChromosome;
  std::vector<Chromosome*> Population;
  std::vector<Chromosome*> NewPopulation;
  std::vector<Cluster*> Clusters;
  double BestFitness;
  double AvgFitness;    // average fitness
  double SigmaFitness;

  std::vector <double> Par;  
  std::vector <std::string> ParName;
  std::vector <double> ParMax;
  std::vector <double> ParMin; 
  std::vector <int> ParCluster;
  std::vector <int> ParaPar;

  SimulatedAnnealing* SA;

  Genetic(void);
  ~Genetic(void);
  void Optimize(bool is_chaining);
  int  GetOptimalInCluster(int no_cluster, double max_fit);
  double ClusterDmax(void);
};

#endif
