//=======================================================================
//            Copyright (c) 2003 by INRA. All rights reserved            
//                Redistribution is not permitted without                
//                the express written permission of INRA                 
//                   eMail : tschiex@toulouse.inra.fr                    
//------------------------------------------------------------------------
// File        : EuGeneTk/Parametrization/LineSearch/Chromosome.h             
// Description : Class Data and Chromosome of the Genetic algorithm
// Authors     : P.Bardou, S.Foissac, M.J.Cros, A.Moisan, T.Schiex       
//=======================================================================

#ifndef CHROMOSOME_H_INCLUDED
#define CHROMOSOME_H_INCLUDED


#include <vector>
#include <string>

class Data {
 public:
  std::vector <double> P;

  Data(int n);
  double GetDistance(Data* d);
  void PutBarycenter(int nb1, Data *data2, int nb2);
  void Recenter(void);
};


class Chromosome {
 public:
  double RawFitness;    // Raw fitness
  double ScaledFitness;    
  bool IsEvaluated;     // true if chromosome already evaluated
  bool IsNew;           // true if new chromosome, false if copy of old one
  int NoCluster;        // Number of the cluster this element belongs to
  bool IsBestInCluster; // true if best element in cluster, false otherwise.  
			// Used only with elitist selection for each cluster
  Data *data;

  Chromosome(int n);
  ~Chromosome(void);
  void Evaluate(void);
  void Print(std::string msg);
};



class Cluster {
 public:
    Data* CentralData;  /* Pointer on the class center of this cluster */
    int NbChrom; /* number of elements in this cluster */

    Cluster(int n);
    ~Cluster(void);
};


#endif
