#ifndef  BACKP_H_INCLUDED
#define  BACKP_H_INCLUDED

// Une liste doublement chainee et circulaire de point de retour
// arriere pour le plus court chemin avec contrainte de duree minimum.
// Le premier element A ne doit pas etre efface. A.Next est le dernier
// insere (le moins couteux de tous), A.Prev est le plus vieux (le
// plus cher).

#include "Const.h"

const unsigned char SwitchStart = 0x1;
const unsigned char SwitchStop  = 0x2;
const unsigned char SwitchAcc   = 0x4;
const unsigned char SwitchDon   = 0x8;
const unsigned char SwitchAny   = 0xF; 

class BackPoint
  {
  private:
    char State;
    unsigned char SwitchType;
    int StartPos;  
    double Cost;
    double Additional;
    BackPoint *Origin;
    
  public:
    BackPoint ();
    BackPoint  (char state, int pos, double cost);
    ~BackPoint();
    void InsertNew(char state, unsigned char Switch,int pos, double cost,BackPoint *Or);
    void Print();
    void Dump();
    void BackTrace(char *Choix);
    void Update(double cost);
    BackPoint *BestUsable(int pos, unsigned char mask, int len,REAL *cost);
    BackPoint *StrictBestUsable(int pos, int len,REAL *cost);
    void Zap();
    BackPoint *Next;
    BackPoint *Prev;    
  };


#endif
