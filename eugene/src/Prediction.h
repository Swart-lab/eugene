// ------------------------------------------------------------------
// Copyright (C) 2004 INRA <eugene@ossau.toulouse.inra.fr>
//
// This program is open source; you can redistribute it and/or modify
// it under the terms of the Artistic License (see LICENSE file).
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
//
// You should have received a copy of Artistic License along with
// this program; if not, please see http://www.opensource.org
//
// $Id$
// ------------------------------------------------------------------
// File:     Prediction.h
// Contents: class Prediction
// ------------------------------------------------------------------

#ifndef  PREDICTION_H_INCLUDED
#define  PREDICTION_H_INCLUDED

#include <cstdio>
#include <vector>
#include <algorithm>

#include "Const.h"
extern "C"{
#include "GDIF/gdIF.h"
}
#include "SensorIF.h"
#include "Param.h"
#include "Hits.h"

//SEB
#include "Prediction_cte.h"
#include "Gff3Line.h"
#include "System.h"
#include <fstream>
#include <iostream>
//SEB

class MasterSensor;
class DNASeq;


/*************************************************************
 **                       Feature Object                    **
 *************************************************************/
class Feature
{
 private:

 public:
  int  number;
  signed char state;
  int  start;
  int  end;
  char strand;
  int  phase;
  int  framegff;
  int  frame;
  
  Feature  ();
  Feature  (signed char, int, int, char, int);
  ~Feature ();
};


/*************************************************************
 **                         Gene Object                     **
 *************************************************************/
class Gene
{
  friend class Prediction;
  
 private:
  int complete;  //+1:init|sngl +2:term|sngl 3:ALL 
  int exNumber;
  int inNumber, inLength;
  int utrLength, mrnaLength, geneLength;
  void Update   ();
  void clear	();
 public:
  int cdsStart, cdsEnd, trStart, trEnd;
  int tuStart, tuEnd; // transcript unit (required for gff3 + alternatives variants)
  int exLength;
  int geneNumber;
  bool isvariant; // true = is a splice variant
  int hasvariant; // has n splice variant (or is number n splice variant if is variant is true)
  std::vector <Feature*> vFea;
  
  Gene  ();
  ~Gene ();

  int isDifferent (const Gene& o, int threshold);
  bool operator== (const Gene& o);
  inline int nbFea() {return vFea.size();};
  void AddFeature(signed char state, int start, int end);
  void PrintInfo (FILE*, int, char*);
  char GetVariantCode(void);

};


/*************************************************************
 **                        Prediction                       **
 ************************************************************/
class Prediction
{
 private:
  MasterSensor *MS;
  DNASeq       *X;
  unsigned char *ESTMatch;

  inline unsigned char GetESTMatch(int i) {return (ESTMatch ? ESTMatch[i] : 0);};
  void  clear			(); 
  void  PrintGff        (FILE*, char*);
//SEB
  void  PrintGff3        (std::ofstream&, char*, char);
//SEB
  void  PrintEgnL       (FILE*, char*, int a=0);
  void  PrintEgnS       (FILE*);
  void  PrintEgnD       (FILE*);
  void  PrintHtml       (FILE*, char*);

  // Convert state in string
  char* State2EGNString (int);
  char* State2GFFString (int);
  
  // Verif coherence EST: calcul le nombre de nuc. coherents et
  // incoherents avec les match est
  // debut/fin/etat: debut et fin de la seq. dont l'etat est etat
  // cons/incons: retour des valeurs
  void  CheckConsistency(int debut, int fin, int etat,
			 int* cons, int* incons);

  // Check and Trim UTR according to EST evidence
  void ESTScan();
  bool UTRCheckAndTrim(int* debut, int* fin, int etat);
  void UpdateAndDelete();

  // PrintHtml : -ph print the start of the HTML output
  void  StartHTML       (char*, FILE*);
  
  // PrintHtml : -ph print the end of the HTML output
  void  EndHTML         (FILE*);

  //SEB
  Gff3Line* fillGff3Line(int type_sofa, int start, int end,
                            char strand, int framegff);
  void setGff3Attributes(Gff3Line* line, int type_egn,
                          int type_sofa, std::string fea_name,
                          int j, char code, std::string gene_id, bool coding);
  bool previousExonMustBeUpdated(Gff3Line* line, int start);
  //SEB


 public:
  std::vector <Gene*> vGene;
  int    nbGene;
  char   seqName[FILENAME_MAX+1];
  double optimalPath;

  Prediction  ();
  Prediction  (int From, int To, std::vector <int> vPos,
         std::vector <signed char> vState);
  ~Prediction ();
  bool IsOriginal(Prediction* optPred, std::vector <Prediction*>& altPreds, int seuil);
  void  TrimAndUpdate (DNASeq*);
  void  DeleteOutOfRange(int s,int e);
  void  Print         (DNASeq*, MasterSensor*, FILE *OTP_OUT=NULL, char append = 0);
  void  PrintGeneInfo (FILE*);
  void  PlotPred      ();
  char  GetStateForPos(int);
  Gene *FindGene(int start, int end);

  // Need by Sensor Tester
  char* IsStart       (int);
  char* IsStop        (int);
  char* IsDon         (int);
  char* IsAcc         (int);
  bool  IsState       (DATA::SigType sig_type, int pos, char strand);
};


// Sensor Contents used by each track
const enum DATA::ContentsType SensorContents[NbTracks] = {
  DATA::ExonF1, DATA::ExonF2, DATA::ExonF3,
  DATA::ExonR1, DATA::ExonR2, DATA::ExonR3,
  DATA::ExonF1, DATA::ExonF2, DATA::ExonF3,
  DATA::ExonR1, DATA::ExonR2, DATA::ExonR3,
  DATA::ExonF1, DATA::ExonF2, DATA::ExonF3,
  DATA::ExonR1, DATA::ExonR2, DATA::ExonR3,
  DATA::ExonF1, DATA::ExonF2, DATA::ExonF3,
  DATA::ExonR1, DATA::ExonR2, DATA::ExonR3,
  DATA::IntronF, DATA::IntronF, DATA::IntronF,
  DATA::IntronR, DATA::IntronR, DATA::IntronR,
  DATA::IntronF, DATA::IntronF, DATA::IntronF,
  DATA::IntronR, DATA::IntronR, DATA::IntronR,
  DATA::InterG, 
  DATA::UTR5F, DATA::UTR3F, 
  DATA::UTR5R, DATA::UTR3R, 
  DATA::IntronUTRF, DATA::IntronUTRR, 
  DATA::IntronUTRF, DATA::IntronUTRR
};

inline int PhaseAdapt(char p) {return State2Frame[(int)p];}


#endif
