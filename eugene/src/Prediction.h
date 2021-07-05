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
// $Id: Prediction.h,v 1.49 2012-11-07 09:44:47 sallet Exp $
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

#include "Prediction_cte.h"
#include "Gff3Line.h"
#include "System.h"
#include "State.h"
#include <fstream>
#include <iostream>

class MasterSensor;
class DNASeq;


/*************************************************************
 **                       Feature Object                    **
 *************************************************************/
class Feature
{
 private:
 public:
  int    number;
  int    start;
  int    end;
  char   strand;
  int    phase;
  int    frame;
  State* featureState;
  int    cds_start; // if feature is a cds, cds start position (can be different from start if there is a FS)
  int    cds_end;   // if feature is a cds, cds end position (can be different from end if there is a FS)
  
  Feature  ();
  Feature  (signed char, int, int);
  ~Feature ();
  bool      Overlap (const Feature&);
  bool      IsCodingExon();
  bool      IsIntron();
  bool      IsUTR();
  bool      IsUTRIntron();
  bool      IsIntergenic();
  bool      IsTranscribedAndUnspliced();
  bool      IsNpcRna();
  bool      IsBicoding();
  void      ComputePhase(int);
  void      SetState(char);
};


/*************************************************************
 **                         Gene Object                     **
 *************************************************************/
class Gene
{
  friend class Prediction;
  
 private:
  int  complete;  //+1:init|sngl +2:term|sngl 3:ALL 
  int  exNumber;
  int  inNumber, inLength;
  int  utrLength, mrnaLength, geneLength;
  char strand;
  bool isNpcRna;   // true if it is a non protein coding rna
  int  operonNb; // operon number
  //std::vector<int> vFrameShift; TODO: save the positions of the frameshift in a array
  int hasFrameShift;
  void Update   (int);
  void clear    ();
 public:
  int cdsStart, cdsEnd, trStart, trEnd;
  int tuStart, tuEnd; // transcript unit (required for gff3 + alternatives variants)
  int exLength;
  int geneNumber;
  bool isvariant; // true = is a splice variant
  int hasvariant; // has n splice variant (or is number n splice variant if is variant is true)
  std::vector <Feature*> vFea;
  
  Gene  ();
  Gene (std::string line, DNASeq& seq);
  ~Gene ();

  int isDifferent (const Gene& o, int threshold);
  bool operator== (const Gene& o);
  bool HasSameExons(const Gene& o);
  bool Overlap (const Gene& o);
  float GetOverlapWith ( const Gene& g) ;
  inline int nbFea() {return vFea.size();};
  int  GetExonNumber();
  char GetStrand();
  int  GetOperonNb();
  void SetOperonNb(int);
  void AddFeature(signed char state, int start, int end);
  void PrintInfo (FILE*, int, char*);
  void Print();
  std::string GetVariantCode(void) const;
  bool IsNpcRna();
  void ManageFrameShift(int);
  void ConvertToNpcRNA(int);

};

/*************************************************************
 **                         Operon Object                     **
 *************************************************************/
class Operon
{
 private:
  int  start;
  int  end;
  char strand;
  int  number;
 public:
  std::vector <Gene*> vGenes;
  Operon();
  Operon(int nb);
  ~Operon();
  int  GetStart();
  int  GetEnd();
  int  GetNumber();
  char GetStrand();
  void SetNumber(int );
  void AddGene(Gene* gene);
  void Print();
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
  void  clear     (); 
  void  PrintGff        (FILE*, char*);
  void  PrintGff3        (std::ofstream&, char*, char);
  void  PrintEgnL       (FILE*, char*, int a=0);
  void  PrintEgnS       (FILE*);
  void  PrintEgnD       (FILE*);
  void  PrintHtml       (FILE*, char*);
  
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

  Gff3Line* fillGff3Line(int type_sofa, int start, int end,
                            char strand, int phasegff);
 void setGff3Attributes(Gff3Line* line, State* featState,
                          int type_sofa, std::string fea_name,
                          int j, std::string code_variant, std::string gene_id, bool coding);
  bool previousExonMustBeUpdated(Gff3Line* line, int start);


 public:
  std::vector <Gene*>   vGene;
  std::vector <Operon*> vOperon;
  int    nbGene;        // Number of predicted genes (including gene non coding for protein)
  char   seqName[FILENAME_MAX+1];
  double optimalPath;

  Prediction  ();
  Prediction  (int From, int To, std::vector <int> vPos,
         std::vector <signed char> vState);
  Prediction  (const std::string & desc, DNASeq*);
  Prediction  (char [FILENAME_MAX+1], DNASeq *);
  ~Prediction ();
  
  void LoadGene(std::vector<GeneFeature*>&, std::vector<int>&, std::vector<signed char>&, DNASeq*); 
  void Init (int From, int To, const std::vector <int>& , const std::vector <signed char>&);
  bool IsOriginal(Prediction* optPred, std::vector <Prediction*>& altPreds, int seuil);
  void  TrimAndUpdate (DNASeq*);
  void  SanityCheck();
  void  DeleteOutOfRange(int s,int e, char strand);
  void  Print         (DNASeq*, MasterSensor*, FILE *OTP_OUT=NULL, char append = 0, char variant = 0);
  void  PrintGeneInfo (FILE*);
  void  PlotPred      ();
  State* GetStateAtPos(int); 
  Gene *FindGene(int start, int end, char strand = 0);
  void Print();
  std::vector<int> Eval(Prediction* ref, int offset, bool onlyCodingGene);
  std::vector<int> EvalGene(Prediction* ref, int start, int end, bool onlyCodingGene);
  std::vector<int> EvalExon(Prediction* ref, int start, int end);
  std::vector<Feature*> GetExons(int start, int end);
  std::vector<Gene*>    GetGenes(int start, int end, bool onlyCodingGene);
  void  ComputeOperons();
  Operon* GetOperon(int nb);
  int GetExonNumber();
  int GetExonLength();
  void AppendPred(const Prediction*);

  // Need by Sensor Tester
  const char* IsStart       (int);
  const char* IsStop        (int);
  const char* IsDon         (int);
  const char* IsAcc         (int);
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
  DATA::IntronF, DATA::IntronF, DATA::IntronF,
  DATA::IntronR, DATA::IntronR, DATA::IntronR,
  DATA::IntronR, DATA::IntronR, DATA::IntronR, DATA::IntronR,
  DATA::InterG, 
  DATA::UTR5F, DATA::UTR3F, 
  DATA::UTR5R, DATA::UTR3R, 
  DATA::IntronUTRF, DATA::IntronUTRR, 
  DATA::IntronUTRF, DATA::IntronUTRR,
  DATA::RNAF, DATA::RNAR,
  DATA::ExonF1, DATA::ExonF1, DATA::ExonF2, // Bicoding: Ces valeurs ne doivent jamais etre appelées
  DATA::ExonR1, DATA::ExonR1, DATA::ExonR2, // Bicoding: Ces valeurs ne doivent jamais etre appelées
  DATA::ExonF1, DATA::ExonF1, DATA::ExonF1, // Bicoding: Ces valeurs ne doivent jamais etre appelées
  DATA::ExonF2, DATA::ExonF2, DATA::ExonF2, // Bicoding: Ces valeurs ne doivent jamais etre appelées
  DATA::ExonF3, DATA::ExonF3, DATA::ExonF3, // Bicoding: Ces valeurs ne doivent jamais etre appelées
  DATA::UIRF, DATA::UIRR
};


#endif
