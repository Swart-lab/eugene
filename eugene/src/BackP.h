#ifndef  BACKP_H_INCLUDED
#define  BACKP_H_INCLUDED

// Une liste doublement chainee et circulaire de point de retour
// arriere pour le plus court chemin avec contrainte de duree minimum.
// Le premier element Path (dans Track) ne doit pas etre
// efface. Path.Next est le dernier insere, Path.Prev est le plus vieux
// ??? A t on vraimen t besoin de Prev ?

#include "Const.h"
#include "Prediction.h"
#include "PenaltyDist.h"

class BackPoint
{
  friend class Track;
  // private:
 public:
  signed char State;
  char Optimal;
  int StartPos;  
  double Cost;
  double Additional;
  BackPoint *Origin;
  
 public:
  
  BackPoint *Next;
  BackPoint *Prev;    
  
  BackPoint ();
  BackPoint  (char state, int pos, double cost);
  ~BackPoint();

  inline void SetState(int i) { State = i; };
  void Print();
};


class Track 
{
 public: 

  int NumBP;
  BackPoint Path;
  PenaltyDist PenD;
  double Optimal;
  int OptPos;

  inline void LoadPenalty(char* name) { PenD.LoadPenaltyDist(name); };
  inline void Update(double cost) {  Path.Next->Additional += cost; Optimal += cost;};
  inline void PayTheSlope() {Path.Next->Additional -= PenD.FinalSlope;
                             Optimal  -= PenD.FinalSlope;  };
  void InsertNew(char state,  int pos, double cost, BackPoint *Or);
  void ForceNew(char state, int pos, double cost, BackPoint *Or);
  BackPoint *BestUsable(int pos, double *cost, int pen = 1);
  Prediction* BackTrace();
  void Dump();
  void Zap();

  Track ();
  ~Track();
  
};

#endif
