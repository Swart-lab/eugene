#ifndef DAG_H_INCLUDED
#define DAG_H_INCLUDED
#endif

#include "Const.h"
#include "Prediction.h"
#include "BackP.h"
#include "SensorIF.h"
#include "Param.h"
#include "DNASeq.h"

class DAG
{
 private: 
  DNASeq *TheSeq;
  int StartPosition;
  int EndPosition;
  double ExPrior, InPrior, IGPrior, FivePrior, ThreePrior, IntronFivePrior;
  int    MinCDSLen;
  double SplicedStopPen;
  int estuse;
  double NormalizingPath;
  Track  LBP[NbTracks];

 public:
  char EvidenceName[FILENAME_MAX+1];//   **
  Prediction *pred; //   **

  DAG();
  DAG (int start, int end, Parameters &PAR, DNASeq *Seq);
  DAG (int start, int end, DAG *RefDag,char* name);
  ~DAG ();
  void WeightThePrior();
  void LoadDistLength();
  void ShortestPathAlgoForward (int position, DATA Data);

  //  void ShortestPathAlgoBackward (int position, DATA Data, int Data_Len, int NoContentsUpdate=0);
  //  int ChooseTheBestTrackIndex (int position, double *bestscore, int sense);
  //  int ChooseTheBestTrackIndex (int position, double *bestscore, DAG* DagR);
  //  void InsertLastBestUsableBP (int position, int sense=1);
  //  void SubstractContents (DATA *data);
  //  double ConnectDag (DAG *RightDag, int position, DATA * data, int datalength);
  double BuildPrediction ();
  void BuildReversePrediction (DAG* dag);
  //  void CleanPrediction (DAG* dag, DAG* dagrev);
  
  inline double NormalizingPathScore() { return NormalizingPath; };
  inline int GetStart() { return StartPosition; };
  inline void NormalizingPathUpdate(double maxi) {NormalizingPath+=maxi;};
  void Print();
};
// **  a passer en private
