//=======================================================================
//            Copyright (c) 2003 by INRA. All rights reserved            
//                Redistribution is not permitted without                
//                the express written permission of INRA                 
//                   eMail : tschiex@toulouse.inra.fr                    
//------------------------------------------------------------------------
// File        : EuGeneTk/Parametrization/LineSearch/Genetic.cc             
// Description : Class of the Genetic algorithm
// Authors     : P.Bardou, S.Foissac, M.J.Cros, A.Moisan, T.Schiex       
//=======================================================================


#include "Genetic.h"

#include "../../EuGene/Param.h"
#include "../ParaOptimization.h"
#include "../Random.h"

extern Parameters PAR;
extern ParaOptimization OPTIM;


//-------------------------------------------------------
// Constructor
//-------------------------------------------------------
Genetic::Genetic (void) : OptiAlgorithm()
{
  int n;
  double d;

  NbRun = PAR.getI("Genetic.NbRun");
  NbGeneration = PAR.getI("Genetic.NbGeneration");
  NbElement = PAR.getI("Genetic.NbElement");

  Rand = new Random(PAR.getC("Genetic.Seed"));

  d = PAR.getD("Genetic.Elitism");
  if (d==0)
    	{cerr <<"SORRY: the old software was just writting no parallel version available when no elitism was asked. This software runs but the functionning has to be valided."<<endl; exit(100);}
  else {
    IsElitist = true;
    Elitism = d;
  }

  n = PAR.getI("Genetic.ScalingType");
  if (n==0)
    IsScaling = false;
  else 
    IsScaling = true;
  Scale = new Scaling (n);

  d = PAR.getD("Genetic.Sharing");
  if (d==0)
    	{cerr <<"SORRY: the old software was just writting no parallel version available when no sharing was asked. This software runs but the functionning has to be valided."<<endl; exit(100);}
  else {
    IsSharing = true;
    Share = new Sharing (d);
  }
  IsClustering =  ( ((string) PAR.getC("Genetic.Clustering")) == "TRUE") ? true : false;

  IsSAMutating =  ( ((string) PAR.getC("Genetic.SA.Mutation")) == "TRUE") ? true : false;
  IsSACrossingOver =  ( ((string) PAR.getC("Genetic.SA.CrossOver")) == "TRUE") ? true : false;

  d_mean=0.0;
  minmaxfactor= 2.0;

  NbCluster = 0;

  // Update Par attributes (the same than Para with just the first elt of clusters IDENTICAL
  int i, q, l;
  vector <bool> FirstEltCluster;
  bool is_to_consider;

  for (i=0; i<ParaClusters.size(); i++) FirstEltCluster.push_back( true );

  for (i=0; i<ParaName.size(); i++) {
    Para.push_back( 0 );

    is_to_consider = false;
    if (ParaClusterRelations[ParaCluster[i]] == IDENTICAL) {
      if ( FirstEltCluster[ParaCluster[i]] ) {
	is_to_consider = true;
	FirstEltCluster[ParaCluster[i]] = false;
      } else {
	// just update the index of the first elt of the cluster which the parameter belongs
	l = 0;
	for (q=0; q<ParCluster.size(); q++) 
	  if (ParaCluster[i] == ParCluster[q]) 
	    {l = q; q = ParCluster.size();}
	ParaPar.push_back( l );
      }
    } else 
	is_to_consider = true;

    if (is_to_consider) {
      Par.push_back( Para[i] );
      ParName.push_back( ParaName[i] );
      ParMin.push_back( ParaMin[i] );
      ParMax.push_back( ParaMax[i] );
      ParCluster.push_back( ParaCluster[i] );

      ParaPar.push_back( Par.size()-1 );
    }
  }

}


//-------------------------------------------------------
// Destructor
//-------------------------------------------------------
Genetic::~Genetic (void)
{
  unsigned int i;

  delete Rand;
  delete CrossOver;
  delete Mutate;
  delete Select;
  delete Scale;
  if (IsSharing) delete Share;
  if (IsSAMutating || IsSACrossingOver) delete SA;

  for (i=0; i<Population.size(); i++) delete Population[i];
  for (i=0; i<NewPopulation.size(); i++) delete NewPopulation[i];
  for (i=0; i<Clusters.size(); i++) delete Clusters[i];
}




//-------------------------------------------------------
// genetic algorithm
//-------------------------------------------------------
void Genetic::Optimize(bool is_chaining)
{
  double max = 0;
  unsigned int i; int j, k;
  int run, gen;

  if (is_chaining)
    {cerr <<"ERROR: ask to chain Genetic Algorithm in Genetic::Optimize."<<endl; exit(100);}


  // not in the constructor because OPTIM.Algorithm[OPTIM.AlgoIndex] is not updated
  CrossOver = new CrossingOver (PAR.getD("Genetic.CrossOverProbability"));
  Mutate = new Mutating (PAR.getD("Genetic.MutationProbability"));
  Select = new Selecting (PAR.getI("Genetic.SelectionType"));
  if (IsSAMutating || IsSACrossingOver)  SA = new SimulatedAnnealing();

 // Display parametrization of the algorithm
  cout <<endl;
  cout << "---------------------------------------------------------------"<<endl;
  cout<< "Optimization of EuGène parameters with the Genetic algorithm"<<endl<<endl;
  cout << "---------------------------------------------------------------"<<endl;
  cout << "Parametrisation of the algorithm:"<<endl<<endl;
  cout << "NbRun: " << NbRun <<endl;
  cout << "NbGeneration: " << NbGeneration <<endl;
  cout << "NbElement: " << NbElement <<endl;
  cout << "CrossOverProbability: " << CrossOver->Proba <<endl;
  cout << "MutationProbability: " << Mutate->Proba <<endl;
  cout << "SelectionType: " << Select->Type <<endl;
  if (!IsElitist) cout <<"No Elitism"; else cout <<"Elitism: " << Elitism; cout <<endl;
  if (!IsSharing) cout << "No sharing"; else cout << "Sharing:" << Share->Value; cout <<endl;
  if (!IsScaling) cout << "No scaling"; else cout << "ScalingType:" << Scale->Type; cout <<endl;
  if (!IsClustering) cout << "No "; cout << "Clustering" <<endl;
  if (!IsSAMutating) cout << "No "; cout << "SA Mutating" <<endl;
  if (!IsSACrossingOver) cout << "No "; cout << "SA CrossingOver" <<endl;
  cout << "Seed: " << Rand->Seed <<endl;
  cout << "Trace: " << ((IsTracing) ? "TRUE" : "FALSE") <<endl<<endl;

  cout << "NbParameter: " << ParaName.size() <<endl;
  cout << "Param: \t"; for (i=0; i<ParaName.size(); i++) cout << ReduceName(ParaName[i]) << "\t"; cout <<endl;
  cout << "Min: \t"; for (i=0; i<ParaMin.size(); i++) cout << ParaMin[i] << "\t"; cout <<endl;
  cout << "Max: \t"; for (i=0; i<ParaMax.size(); i++) cout << ParaMax[i] << "\t"; cout <<endl<<endl<<endl;


  for (k=0; k<NbElement; k++)
    Clusters.push_back( new Cluster(ParName.size()) );
 

  // general loop
  for (run=0; run<NbRun; run++) {
    InitPopulation();

    for (gen=0; gen<NbGeneration; gen++) {
      cout <<endl<< "Run: " <<run<< "\tGen: " << gen 
	   << "--------------------------------------------"<<endl;
      EvalPopulation() ;

      if ((IsSharing) && (IsClustering)) NbCluster = MakeClusters();

      Scale->Scalestat(gen);

      if ((gen==0) && ((IsSACrossingOver)||(IsSAMutating)))   SA->SAInit();

      if (IsSharing) Share->Share() ;

     if (BestFitness == 0.0)
	cout << "No admissible element." <<endl;
      else {
	cout << "Fitness: " << BestChromosome->RawFitness 
	     << "\tAvgFitness: " << AvgFitness
	     << "\tSigmaFitness: " << SigmaFitness <<endl;
	for (i=0; i<BestChromosome->data->P.size(); i++) 
	  cout << BestChromosome->data->P[i] << "\t"; 
	cout<<endl;
	if (BestChromosome->RawFitness > 999999999999.0) gen=NbGeneration;
      }

      /* Gather every element which are best of cluster and mark */
      /* them with the flag best of cluster */
      if ((IsSharing) && (IsClustering))   MarkBestOfClusters();
 
      if ((IsSharing) && (IsClustering) && (gen==(NbGeneration-1))) {
	max = BestChromosome->RawFitness * Share->Value ;
	/* Get elements of the population that are the best of */
	/* each optimal cluster, then print optimal clusters   */
	/* and corresponding best fitness */
	cout <<endl<< "Optimum for the run "<< run <<"-------------------------------------"<<endl;
	for (j=0; j<NbCluster; j++) {
	  /* Put index of best element from cluster j in k if this */
	  /* cluster is optimal else put -1 in it */
	  if (Clusters[j]->NbChrom!=0) {
	    k = GetOptimalInCluster(j, max) ;
	    if (k != -1) {
	      cout << "Cluster " << j << ": ";
	      for (i=0; i<Population[k]->data->P.size(); i++) 
		cout << Population[k]->data->P[i] << "\t"; 
	      cout << "Fitness="<<Population[k]->RawFitness <<"\t";
	      if (BestChromosome == Population[k]) {
		cout <<"(Best element)";
		// Remenber the optimum of the run if it is better than the previous ones
		if ((ParOpti.size()==0) || (BestFitness > FitnessOpti))
		  { ParOpti = Population[k]->data->P; FitnessOpti = Population[k]->RawFitness;}
	      }
	      cout <<endl;
	    }
	  }
	}
      }
    } // gen

    /* If last generation then everything is finished so skip following code */
    if (gen < NbGeneration-1) {
      if ((IsSACrossingOver)||(IsSAMutating)) 	SA->SAUpdate(gen);
	
      // For next generation, we select elements then crossover and mutate them
      Select->Selection();

      CrossOver->Crosseval();

      Mutate->Muteval();
    } 
    cout <<endl;
  } /*run*/

  // Put the best optimum of the runs in Para
  Par = ParOpti; BestFitness = FitnessOpti;
  for (int q=0; q<Para.size(); q++) Para[q] = Par[ParaPar[q]];
  if (run>0) {
    cout <<endl<< "---------------------------------------------------------------"<<endl;
    cout <<"Best optimum for the runs"<<endl;
    for (i=0; i<ParaName.size(); i++) cout << ReduceName(ParaName[i]) << "\t"; cout <<endl;
    for (i=0; i<ParaName.size(); i++) cout << Para[i] << "\t"; 
    cout << "Fitness=" << BestFitness <<endl;
    cout << "---------------------------------------------------------------"<<endl;
  }
}



//-------------------------------------------------------
// Init the Population of chromosome
//-------------------------------------------------------
void Genetic::InitPopulation (void)
{
  Chromosome* c;

  for (int i=0; i<NbElement; i++) {
    c = new Chromosome(ParName.size());
    Population.push_back(c);
    c = new Chromosome(ParName.size());
    NewPopulation.push_back(c);
  }
}


//-------------------------------------------------------
// Eval the Population of chromosome
//-------------------------------------------------------
void Genetic::EvalPopulation (void)
{
  double maxfit = -DBL_MAX;

  for (int i=0; i<NbElement; i++) {
    if (!Population[i]->IsEvaluated) Population[i]->Evaluate();

    /* Update BestChromosome */
    if (Population[i]->RawFitness > maxfit) {
      BestChromosome = Population[i];
      maxfit = Population[i]->RawFitness ;
    }
  }
}


//-------------------------------------------------------
// Update clusters and return number of effective clusters
//-------------------------------------------------------
int Genetic::MakeClusters (void)
{
  int i,j,cl;
  double d,distance, dmin, dmax, sumd, nb_moy, max;
  int nopt, nb_dist, nb_clusters;
  bool cont;

  dmax = d_mean/minmaxfactor;
  dmin=dmax/3.0;
  nb_clusters=0;
  sumd=0;
  nb_dist=0;
  
  for (i=0; i<NbElement; i++) {
    distance=DBL_MAX;
    for (j=0; j<nb_clusters; j++) {
      /* 0 means that this cluster has been merged with an other cluster
	 and is no longer valid*/
      if (Clusters[j]->NbChrom!=0) {
	d = Clusters[j]->CentralData->GetDistance( Population[i]->data );
	sumd+=d; nb_dist++;
	if (d<distance)	{
	  cl=j;
	  distance=d;
	}
      }
    }
     
    if (distance<dmax) {
      Clusters[cl]->CentralData->PutBarycenter(Clusters[cl]->NbChrom,
					Population[i]->data, 1);
      Clusters[cl]->NbChrom++;
      do {
	cont = false;
	for (j=0; j<nb_clusters; j++) {
	  if ( (Clusters[j]->NbChrom!=0) && (j!=cl) ) {
	    d = Clusters[j]->CentralData->GetDistance( Clusters[cl]->CentralData );
	    if (d<dmin)	{
	      cl = MergeClusters(j,cl,i);
	      cont = true;
	      break;
	    }
	  }
	}
      } while (cont);
    } else { 
      Clusters[nb_clusters]->CentralData->P = Population[i]->data->P;
      Clusters[nb_clusters]->NbChrom = 1;
      cl=nb_clusters;
      nb_clusters++;
      IndexesOfBestInClusters.push_back(i);
    }
    Population[i]->NoCluster=cl;
  }

  nb_moy=0.0;
  /* Count number of optimal clusters */
  max = BestChromosome->RawFitness * Share->Value;  
  nopt = CountNopt(nb_clusters, max) ;

  cout << "Nb of  cluster: " << nb_clusters << "\t";
  cout << "Number of optimal clusters: " << nopt <<endl;
  cout << "Nb of element by cluster: ";
  for (i=0; i<nb_clusters; i++) {
    nb_moy += Clusters[i]->NbChrom;
    cout << Clusters[i]->NbChrom <<" ";
    }
  cout <<endl;



  cout << "Mean number of element by clusters: " << nb_moy/nb_clusters << "\t";
  d_mean=sumd/nb_dist;
  cout << "Mean distance:" << d_mean <<"\t";

  max=(double)(nopt)/(double)(nb_clusters) ;
  if ((max>0.85)&&(minmaxfactor<100.0)) {
    minmaxfactor *= 1.05;
    cout << "Higher minmaxfactor= " << minmaxfactor <<endl;
  }
  if ((max<0.75)&&(minmaxfactor>1.0)) {    
    minmaxfactor *= 0.95;
    cout << "Lower  minmaxfactor= " << minmaxfactor <<endl;
  }

  return nb_clusters;
}


//-------------------------------------------------------
// Put best_of_cluster flag to every element which is best within its cluster
//-------------------------------------------------------
void Genetic::MarkBestOfClusters(void) 
{
  int i;
  double Max[NbCluster];

  for (i=0; i<NbCluster; i++)
    Max[i] = 0;

  for (i=0; i<NbElement; i++)
    if (Population[i]->RawFitness > Max[Population[i]->NoCluster]) {
      Max[Population[i]->NoCluster] = Population[i]->RawFitness ;
      IndexesOfBestInClusters[Population[i]->NoCluster] = i ;
    }
  
  for (i=0; i<NbCluster; i++)
    Population[IndexesOfBestInClusters[i]]->IsBestInCluster = true;
}


//-------------------------------------------------------
/* Return index of best element from cluster 'cluster' */
/* if this cluster is optimal i.e the fitness of this  */
// element is greater than max. 
// Return -1 if cluster is not optimal 
//-------------------------------------------------------
int Genetic::GetOptimalInCluster(int no_cluster, double max_fit)
{
  if (Population[IndexesOfBestInClusters[no_cluster]]->RawFitness > max_fit)
    return (IndexesOfBestInClusters[no_cluster]) ;
  else
    return (-1) ;
}


//-------------------------------------------------------
/* ClusterMake divides the population into clusters. The principle of */
/* the algorithm is simple. Each population element is compared to */
/* the existing cluster centers. It is then put into the closest */
/* cluster if its distance to this cluster is less than dmax. If its */
/* distance to the closest cluster is larger than dmax, then a new */
/* cluster is created, and the current population element becomes the */
/* center of the new cluster. If, when adding a population element to */
/* a cluster, the distance of this cluster to another cluster becomes */
/* less than dmin, then the two clusters are merged. */
//-------------------------------------------------------
double Genetic::ClusterDmax(void)
{
  return d_mean/minmaxfactor;
}



//-------------------------------------------------------
/* MergeClusters takes two cluster numbers and merges them into a single cluster.
   The new cluster is the barycenter of the two late clusters.
   One of the cluster is marked as invalid (the number of elements is set to 0).
   The population is updated in order to take into account the disappearance of
   one of the clusters. */
//-------------------------------------------------------
int Genetic::MergeClusters(int cl1, int cl2,int last)
{
  int i;
    
  if (Clusters[cl1]->NbChrom < Clusters[cl2]->NbChrom) {
    i=cl1; cl1=cl2; cl2=i;
  }
  Clusters[cl1]->CentralData->PutBarycenter(Clusters[cl1]->NbChrom,
				     Clusters[cl2]->CentralData, Clusters[cl2]->NbChrom);
  Clusters[cl1]->NbChrom += Clusters[cl2]->NbChrom;
  Clusters[cl2]->NbChrom = 0;
  
  for (i=0; i<last; i++) {
    if (Population[i]->NoCluster == cl2)
      Population[i]->NoCluster = cl1;
  }
  return cl1;
}


//-------------------------------------------------------
//-------------------------------------------------------
int Genetic::CountNopt(int nb_clusters, double max) 
{
  int i, nopt=0 ;
  int optimal_cluster[nb_clusters] ;

  /* At the beginning none of the clusters is optimal */
  for (i=0; i<nb_clusters; i++)
    optimal_cluster[i]=0 ;

  for (i=0; i<NbElement; i++)
    if (Population[i]->RawFitness > max)
      optimal_cluster[Population[i]->NoCluster]=1 ;

  /* Collect final results */
  for (i=0; i<nb_clusters; i++)
    nopt += optimal_cluster[i];

  return (nopt) ;
}


