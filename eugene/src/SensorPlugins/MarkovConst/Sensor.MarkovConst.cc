#include "Sensor.MarkovConst.h"

/*************************************************************
 **                      SensorMarkovConst                   **
 *************************************************************/

extern Parameters PAR;

// ----------------------
//  Default constructor.
// ----------------------
SensorMarkovConst :: SensorMarkovConst (int n, DNASeq *X) : Sensor(n)
{
  type = Type_Content;
}

// ----------------------
//  Default destructor.
// ----------------------
SensorMarkovConst :: ~SensorMarkovConst ()
{
}

// ----------------------
//  Init Markov.
// ----------------------
void SensorMarkovConst :: Init (DNASeq *X)
{
  transCodant = PAR.getD("MarkovConst.Coding*"); //Exon
  transIntron = PAR.getD("MarkovConst.Intron*"); //IntronF
  transInter  = PAR.getD("MarkovConst.Inter*");  //InterG
  transUTR5   = PAR.getD("MarkovConst.UTR5*");   //UTR5
  transUTR3   = PAR.getD("MarkovConst.UTR3*");   //UTR3

  minGC = PAR.getD("MarkovConst.minGC",GetNumber());
  maxGC = PAR.getD("MarkovConst.maxGC",GetNumber());
}

// -----------------------
//  GiveInfo Content Markov.
// -----------------------
void SensorMarkovConst :: GiveInfo (DNASeq *X, int pos, DATA *d)
{
  if ((X->Markov0[BitG] + X->Markov0[BitC]) > minGC &&
      (X->Markov0[BitG] + X->Markov0[BitC]) <= maxGC)
    {
      for(int i=0;i<6;i++)
	d->contents[i] += log(transCodant); //Exon
      
      d->contents[6] += log(transIntron); //IntronF
      d->contents[7] += log(transIntron); //IntronR
      d->contents[8] += log(transInter);  //InterG
      d->contents[9] += log(transUTR5);   //UTR5'F
      d->contents[10]+= log(transUTR5);   //UTR5'R
      d->contents[11]+= log(transUTR3);   //UTR3'F
      d->contents[12]+= log(transUTR3);   //UTR3'R
    }
}

void AmplifyScore(double Score [], unsigned int normopt)
{
  double Sum = 0.0;
  double Min = -log(0.0);
  int i;
  
  for  (i = 0;  i < 9;  i ++) 
    if (Min > Score[i]) Min = Score[i];
  
  switch (normopt) {
  case 0:
    for  (i = 0;  i < 9;  i ++) 
      Score[i] = exp(Min-Score[i]);
    // pas de normalisation
    
    break;
    
  case 1:
    // on suppose que une phase au plus peut coder
    for  (i = 0;  i < 9;  i ++) {
      Score[i] = exp(Score[i]-Min);
      Sum += Score[i];
    }
    
    for  (i = 0;  i < 9;  i ++) 
      Score [i] /= Sum;
    break;
    
  case 2:
    // chaque phase peut coder independamment
    for (i = 0; i < 9; i++) Score[i] = exp(Score[i]-Min);
    for (i = 0; i < 6; i++) {
      Sum += Score[i];
      Score [i] /= (Score[6]+Score[7]+Score[8]+Score[i]);
   }
    
    Score[6] /= (Sum+Score[6]+Score[7]+Score[8]);
    Score[7] /= (Sum+Score[6]+Score[7]+Score[8]);
    Score[8] /= (Sum+Score[6]+Score[7]+Score[8]);
    break;
   }
  return;
}

// ----------------------------
//  Plot Sensor information
// ----------------------------
void SensorMarkovConst :: Plot(DNASeq *TheSeq)
{
}

// ------------------
//  Post analyse
// ------------------
void SensorMarkovConst :: PostAnalyse(Prediction *pred)
{
}
