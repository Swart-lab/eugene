//=======================================================================
//            Copyright (c) 2003 by INRA. All rights reserved            
//                Redistribution is not permitted without                
//                the express written permission of INRA                 
//                   eMail : tschiex@toulouse.inra.fr                    
//------------------------------------------------------------------------
// File        : EuGeneTk/Parametrization/ParaOptimization.cc
// Description : The ParaOptimization class optimizes EuGene parameters
// Authors     : P.Bardou, S.Foissac, M.J.Cros, A.Moisan, T.Schiex       
//=======================================================================

#include <iostream>
#include <fstream>

#include "ParaOptimization.h"

#include "../EuGene/Param.h"
#include "../EuGene/Prediction.h"
#include "../EuGene/DNASeq.h"
#include "../EuGene/MSensor.h"
#include "../EuGene/Output.h"
#include "LineSearch/LineSearch.h"
#include "Genetic/Genetic.h"

extern Parameters PAR;
extern MasterSensor*    MS;

extern Prediction* Predict (DNASeq* TheSeq, MasterSensor* MSensor);


//-------------------------------------------------------
// Destructor
//-------------------------------------------------------
ParaOptimization::~ParaOptimization(void)
{
  unsigned int i;
  // delete created instances 
  for (i=0; i<Algorithms.size(); i++) delete Algorithms[i]; 
  for (i=0; i<Sequences.size(); i++) delete Sequences[i];
  for (i=0; i<MSensors.size(); i++) delete MSensors[i];
}


//-------------------------------------------------------
// ParaOptimize : Parameters optimization
//-------------------------------------------------------
void ParaOptimization::ParaOptimize (int argc, char * argv [])
{ 
  std::string filename;
  OptiAlgorithm* algo;
  bool is_chaining = false;

  // Initialisation
  Init(argc, argv);

  // Optimize parameters
  for (unsigned int i=0; i<Algorithms.size(); i++) {
    AlgoIndex = i;
    if (i>0) is_chaining = true;
    Algorithms[i]->Optimize(is_chaining);
  }

  // Write new parameters file
  algo = Algorithms.back();
  filename = PAR.WriteParam( ExecutableName.c_str(), algo->ParaName, algo->Para);
  std::cerr <<std::endl << "A new parameter file " << filename << " is written." << std::endl;
}


//-------------------------------------------------------
// Initialisation 
//-------------------------------------------------------
void ParaOptimization::Init(int argc, char * argv [])
{
  FILE   *fp;
  std::string algo_name;

  ExecutableName = argv[0];
  IsTest = (((std::string) PAR.getC("ParaOptimization.Test") == "TRUE") ? true : false);
  // Inhibit graphic mode
  PAR.set("Output.graph", "FALSE");

  // Creation of required instances of optimization algorithms  
  algo_name = PAR.getC("ParaOptimization.Algorithm");
  if (algo_name == "GENETIC") 
    Algorithms.push_back( new Genetic() );
  else
    if (algo_name == "LINESEARCH")  
      Algorithms.push_back( new LineSearch() );
    else
      if (algo_name == "GENETIC+LINESEARCH") {
	Algorithms.push_back( new Genetic() );
	Algorithms.push_back( new LineSearch() );
      } else {
	std::cerr <<"ERROR: Bad optimization algorithm "<<algo_name<<" in the parameter file"<<std::endl;
	exit(100);
      }

  if (!IsTest) {
    int sequence;
    TrueCoordFile = PAR.getC("ParaOptimization.TrueCoordFile");

    // Update the sequences list
    std::cout << "Loading sequence(s) file(s) ...";
    for (sequence = optind; sequence < argc ; sequence++) {
      fp = (*argv[sequence] ? FileOpen (NULL, argv[sequence], "r") : stdin);    
      if (fp == NULL) {
	std::cerr <<"ERROR: Cannot open fasta file "<<argv[sequence]<<std::endl; 
	exit(100);
      }        
      Sequences.push_back( new DNASeq(fp) );    
      if (fp != stdin) fclose(fp);
      SeqNames.push_back(argv[sequence]);
    }
    std::cout << "done (" << Sequences.size() << " sequence(s))" << std::endl;
  
    // Update the master sensors list
    for (sequence=0; sequence<(int)Sequences.size(); sequence++) {
      PAR.set("fstname", SeqNames[sequence].c_str());
      MSensors.push_back(new MasterSensor);
      MS = MSensors[sequence];
      MSensors[sequence]->InitMaster(Sequences[sequence]);
    }
  }
}


//-------------------------------------------------------
// ParaEvaluate : Evaluate the parameters Algorithms[AlgoIndex]->Para
//-------------------------------------------------------
double ParaOptimization::ParaEvaluate (bool is_detail_required)
{
  OptiAlgorithm* algo = Algorithms[AlgoIndex];

  double fitness = 0;
  unsigned int i;
  std::string cmde, fic_pred, fic_eval, eval_options;
  std::ifstream feval;
  double spg, sng, spe, sne;
  Prediction* pred;
  char* c = "sequence"; char ** cc = &c;
  FILE   *fp;

  DetailedEvaluation = "";

  if (algo->Para.size() > 0) {
    if (!IsTest) {

      // Update new value for parameters
      for (i=0; i<Algorithms[AlgoIndex]->ParaName.size(); i++) 
	PAR.setD(Algorithms[AlgoIndex]->ParaName[i].c_str(), Algorithms[AlgoIndex]->Para[i]);

      // update sensors
      for (i=0; i<Sequences.size(); i++) {
	PAR.set("fstname", SeqNames[i].c_str());
	MS = MSensors[i]; 
	MSensors[i]->InitSensors(Sequences[i]);
      }

      fic_pred = "tmp%pred%" + ExecutableName;
      cmde = "rm -f " + fic_pred; system(cmde.c_str());
      fp = fopen(fic_pred.c_str(),"w");
      for (i=0; i<Sequences.size(); i++) {
	PAR.set("fstname", SeqNames[i].c_str());
	pred = Predict(Sequences[i], MSensors[i]);  
	Output(Sequences[i], MSensors[i], pred, 1, 0, cc, fp);
	if (i!=Sequences.size()-1) fprintf(fp,"\n");
	delete pred;
      }
      fclose(fp);
      fic_eval = "tmp%eval%" + ExecutableName;
      if (!is_detail_required) eval_options = " -ps -o1";
      else eval_options = " -o1";
      cmde = "rm -f " + fic_eval; system(cmde.c_str());
      cmde =  "../Procedures/Eval/evalpred.pl " + TrueCoordFile + " " +
	 fic_pred + eval_options + " > " + fic_eval; system(cmde.c_str());

      feval.open(fic_eval.c_str(), std::ios::in);
      if (!is_detail_required) feval >> spg >> sng >> spe >> sne;
      else { 
	char line[256]; 
	for (i=0; i<4; i++) feval.getline(line,256); // ignore the 4th lines
	for (i=0; i<3; i++) {feval.getline(line,256);DetailedEvaluation=DetailedEvaluation+"** "+line+"\n";}
      }
      feval.close(); 

      if (!is_detail_required)  fitness = pow(spg*sng*spe*sne,0.25);
      else fitness = 0;
    } else {
      fitness = 1;
      for (i=0; i<algo->Para.size(); i++)
	fitness *= NormalLaw( algo->Para[i] );
    }
  }
  return fitness;
}


//-------------------------------------------------------
// Normal Law 
//-------------------------------------------------------
double ParaOptimization::NormalLaw(double x)
{
  return exp(-pow(x-0.5,2)/2) / sqrt(2*M_PI);
}


