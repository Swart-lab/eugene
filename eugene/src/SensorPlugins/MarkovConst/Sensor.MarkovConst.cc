#include "Sensor.MarkovConst.h"

/*************************************************************
 **                      SensorMarkovConst                   **
 *************************************************************/

extern Parameters PAR;

// ----------------------
//  Default constructor.
// ----------------------
SensorMarkovConst :: SensorMarkovConst (int n) : Sensor(n)
{
  transCodant = PAR.getD("Markov.transCodant"); //Exon
  transIntron = PAR.getD("Markov.transIntron"); //IntronF
  transInter  = PAR.getD("Markov.transInter");  //InterG
  transUTR5   = PAR.getD("Markov.transUTR5");   //UTR5
  transUTR3   = PAR.getD("Markov.transUTR3");   //UTR3

  minGC = PAR.getD("Markov.minGC",GetNumber());
  maxGC = PAR.getD("Markov.maxGC",GetNumber());
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
  type = Type_Content;
  
  if(PAR.getI("Output.graph")) Plot(X);
  value=0.25;
}

// -----------------------
//  GiveInfo Content Markov.
// -----------------------
void SensorMarkovConst :: GiveInfo (DNASeq *X, int pos, DATA *d)
{
  //  for(int i=0;i<13;i++)
  //    d->ContentScore[i] += log(value);

  if ((X->Markov0[BitG] + X->Markov0[BitC]) > minGC &&
      (X->Markov0[BitG] + X->Markov0[BitC]) <= maxGC)
    {
      for(int i=0;i<6;i++)
	d->ContentScore[i] += log(value) + log(transCodant); //Exon
      
      d->ContentScore[6] += log(value) + log(transIntron); //IntronF
      d->ContentScore[7] += log(value) + log(transIntron); //IntronR
      d->ContentScore[8] += log(value) + log(transInter);  //InterG
      d->ContentScore[9] += log(value) + log(transUTR5);   //UTR5'F
      d->ContentScore[10]+= log(value) + log(transUTR5);   //UTR5'R
      d->ContentScore[11]+= log(value) + log(transUTR3);   //UTR3'F
      d->ContentScore[12]+= log(value) + log(transUTR3);   //UTR3'R
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
  int window, normopt;
  DATA data;
  double *Score = data.ContentScore;
  double NScore[9], LScore[9] = {-1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0,-1.0,-1.0};  
  int i,j,p;

  window = PAR.getI("Output.window")*2+1;
  normopt = PAR.getI("Output.normopt");
  
  for (j = 0 ; j < 9 ; j++)
    NScore[j] = Score[j] = 0.0;
  
  for (i = 0; i < window/2; i++) {
    GiveInfo(TheSeq,i, &data);
    for (j = 0 ; j < 9 ; j++) {
      NScore[j] += Score[j];
      Score[j] = 0;
    }
  }

  for (i=0; i<TheSeq->SeqLen; i++)  {

    if (i-window/2 >= 0) {
      GiveInfo(TheSeq,i-window/2, &data);
      for (j = 0 ; j < 9 ; j++) {
	NScore[j] -= Score[j];
	Score[j] = 0;
      }
    }

    if (i+window/2 < TheSeq->SeqLen) {
      GiveInfo(TheSeq,i+window/2,&data);
      for (j = 0 ; j < 9 ; j++) 
	NScore[j] += Score[j];
    }
    
    for (j = 0 ; j < 9 ; j++) Score[j] = NScore[j];
    AmplifyScore(Score,normopt);
    
    if (LScore[0] < 0) for (j = 0 ; j < 9 ; j++) LScore[j] = Score[j];

    p = ((i == 0) ? 0 : i-1);
    
    PlotLine(p,i, 1, 1,LScore[0],Score[0],3);
    PlotLine(p,i, 2, 2,LScore[1],Score[1],3);
    PlotLine(p,i, 3, 3,LScore[2],Score[2],3);
    PlotLine(p,i,-1,-1,LScore[3],Score[3],3);
    PlotLine(p,i,-2,-2,LScore[4],Score[4],3);
    PlotLine(p,i,-3,-3,LScore[5],Score[5],3);
    PlotLine(p,i, 4, 4,LScore[6],Score[6],3);
    PlotLine(p,i,-4,-4,LScore[7],Score[7],3);
    PlotLine(p,i, 0, 0,LScore[8],Score[8],3);
      
    for (j = 0 ; j < 9 ; j++) {
      LScore[j] = Score[j];
      Score[j] = 0;
    }
  }
}

// ------------------
//  Post analyse
// ------------------
void SensorMarkovConst :: PostAnalyse(Prediction *pred)
{
}
