#include "Sensor.MarkovIMM.h"

/*************************************************************
 **                     SensorMarkovIMM                     **
 *************************************************************/
extern Parameters PAR;

// ---------------------------------------------------- 
// Normalisation sur les 6 phases + non codant
// ---------------------------------------------------- 

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

// ----------------------
//  Default constructor.
// ----------------------
SensorMarkovIMM :: SensorMarkovIMM (int n) : Sensor(n)
{
  FILE *fp;
  int i;
  
  transCodant = PAR.getD("EuGene.transCodant"); //Exon
  transIntron = PAR.getD("EuGene.transIntron"); //IntronF
  transInter  = PAR.getD("EuGene.transInter");  //InterG
  transUTR5   = PAR.getD("EuGene.transUTR5");   //UTR5
  transUTR3   = PAR.getD("EuGene.transUTR3");   //UTR3

  minGC = PAR.getD("MarkovIMM.minGC",GetNumber());
  maxGC = PAR.getD("MarkovIMM.maxGC",GetNumber());
  
  if (! (fp = FileOpen(EugDir, PAR.getC("MarkovIMM.matname",GetNumber()), "rb"))) {
    fprintf(stderr, "cannot open matrix file %s\n", PAR.getC("matname"));
    exit(2);
  }
  
  fprintf(stderr,"Loading IMM...");
  fflush(stderr);
  
  // On essaie d'abord de charger les 5 modeles fondamentaux (3 ex/int/interG)
  for  (i = 0;  i < 5;  i ++) {
    IMMatrix[i] = new BString_Array(MODEL_LEN, ALPHABET_SIZE);
    if (IMMatrix[i]->Read(fp)) {
      fprintf(stderr,"Model %d unreadable in %s. Aborting.\n",i+1,PAR.getC("matname"));
      exit(1);
    } 
    fprintf(stderr,"%d ",i+1);
    fflush(stderr);
  }

  // On essaie ensuite de lire un 6eme modele. Si cela echoue,
  // le modele intronique est utilise pour les UTR
  IMMatrix[6] = new BString_Array(MODEL_LEN, ALPHABET_SIZE);
  if (IMMatrix[6]->Read(fp)) {
    fprintf(stderr,"- No UTR model found, using introns model. ");
    delete IMMatrix[6];
    IMMatrix[6] = IMMatrix[3];
    IMMatrix[5] = IMMatrix[3];
  } else {
    fprintf(stderr,"6 ");
    IMMatrix[5] = new BString_Array(MODEL_LEN, ALPHABET_SIZE);
    if (IMMatrix[5]->Read(fp)) {
      fprintf(stderr,"- No second UTR model found, using intron model. ");
      delete IMMatrix[5];
      IMMatrix[5] = IMMatrix[3];
    } else fprintf(stderr,"7 ");
  }  
  fprintf(stderr,"done\n");
  fclose(fp);
}

// ----------------------
//  Default destructor.
// ----------------------
SensorMarkovIMM :: ~SensorMarkovIMM ()
{
  // free remaining used memory
  if (IMMatrix[5] != IMMatrix[3]) delete  IMMatrix[5];
  if (IMMatrix[6] != IMMatrix[3]) delete  IMMatrix[6];
  
  for  (int i = 0;  i < 5;  i ++)
    delete  IMMatrix[i];
}

// ----------------------
//  Init MarkovIMM.
// ----------------------
void SensorMarkovIMM :: Init (DNASeq *X)
{
  type = Type_Content;
  
  if(PAR.getI("Output.graph"))
    Plot(X);
}

// -----------------------
//  ResetIter.
// -----------------------
void SensorMarkovIMM :: ResetIter ()
{
}

// ----------------------------
//  GiveInfo Content MarkovIMM.
// ----------------------------
void SensorMarkovIMM :: GiveInfo(DNASeq *X, int pos, DATA *d)
{
  int  Rev,FModelLen;
  int indexF,indexR;

  // If the model is not in its GC% area, simply do nothing
  if ((X->Markov0[BitG] + X->Markov0[BitC]) <= minGC ||
      (X->Markov0[BitG] + X->Markov0[BitC]) > maxGC)
    return;

  // else if the current nuc. is unknow, again do nothing.
  if ((*X)[pos] == 'n') 
    return;

  // Compute possible maximum order
  Rev = X->SeqLen - pos;
  FModelLen = Min(pos+1,MODEL_LEN);
  
  // and indexes in the IMMatrix models
  indexF = (*IMMatrix[0]).AntiString_To_Sub(X,pos,Min(Rev,MODEL_LEN));
  indexR = (*IMMatrix[0]).String_To_Sub(X,pos-FModelLen+1,FModelLen);

  // Exons F
  d->ContentScore[0] += log((double)(*IMMatrix[2-((pos+2)%3)])[indexF]/65535.0);
  d->ContentScore[1] += log((double)(*IMMatrix[2-((pos+1)%3)])[indexF]/65535.0);
  d->ContentScore[2] += log((double)(*IMMatrix[2-(pos%3)])[indexF]/65535.0);
    
  // Exons R
  d->ContentScore[3] += log((double)(*IMMatrix[2-((Rev+1)%3)])[indexR]/65535.0);
  d->ContentScore[4] += log((double)(*IMMatrix[2-(Rev%3)])[indexR]/65535.0);
  d->ContentScore[5] += log((double)(*IMMatrix[2-((Rev+2)%3)])[indexR]/65535.0);
  
  // Introns F/R
  d->ContentScore[6] += log((double)(*IMMatrix[3])[indexF]/65535.0);
  d->ContentScore[7] += log((double)(*IMMatrix[3])[indexR]/65535.0);
    
  // InterG
  d->ContentScore[8] += 
    log(((double)(*IMMatrix[4])[indexF] + (double)(*IMMatrix[4])[indexR])/131071.0);
    
  // UTR 5' F/R
  d->ContentScore[9] += log((double)(*IMMatrix[6])[indexF]/65535.0);
  d->ContentScore[10] += log((double)(*IMMatrix[6])[indexR]/65535.0);
      
  // UTR 3' F/R
  d->ContentScore[11] += log((double)(*IMMatrix[5])[indexF]/65535.0);
  d->ContentScore[12] += log((double)(*IMMatrix[5])[indexR]/65535.0);
  
  for(int i=0; i<6; i++)
    d->ContentScore[i] += log(transCodant); //Exon
  d->ContentScore[6] += log(transIntron); //IntronF
  d->ContentScore[7] += log(transIntron); //IntronR
  d->ContentScore[8] += log(transInter);  //InterG
  d->ContentScore[9] += log(transUTR5);   //UTR5'F
  d->ContentScore[10] += log(transUTR5);   //UTR5'R
  d->ContentScore[11] += log(transUTR3);   //UTR3'F
  d->ContentScore[12] += log(transUTR3);   //UTR3'R
  return;
}

// ------------------------------
//  GiveInfoAt Content MarkovIMM.
// ------------------------------
void SensorMarkovIMM :: GiveInfoAt (DNASeq *X, int pos, DATA *d)
{
  GiveInfo(X, pos, d);
}

// ----------------------------
//  Plot Sensor information
// ----------------------------
void SensorMarkovIMM :: Plot(DNASeq *TheSeq)
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
      GiveInfoAt(TheSeq,i+window/2,&data);
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
void SensorMarkovIMM :: PostAnalyse(Prediction *pred)
{
}
