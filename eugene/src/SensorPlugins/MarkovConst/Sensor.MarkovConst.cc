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
  transCodant  = PAR.getD("MarkovConst.Coding*"); //Exon
  transIntron  = PAR.getD("MarkovConst.Intron*"); //Intron
  transIntronU = PAR.getD("MarkovConst.IntronUTR*"); //IntronUTR
  transInter   = PAR.getD("MarkovConst.Inter*");  //InterG
  transUTR5    = PAR.getD("MarkovConst.UTR5*");   //UTR5
  transUTR3    = PAR.getD("MarkovConst.UTR3*");   //UTR3

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
      
      d->contents[DATA::IntronF] += log(transIntron);
      d->contents[DATA::IntronR] += log(transIntron);
      d->contents[DATA::InterG] += log(transInter);
      d->contents[DATA::UTR5F] += log(transUTR5);
      d->contents[DATA::UTR5R]+= log(transUTR5);
      d->contents[DATA::UTR3F]+= log(transUTR3);
      d->contents[DATA::UTR3R]+= log(transUTR3);
      d->contents[DATA::IntronUTRF] += log(transIntronU);
      d->contents[DATA::IntronUTRR] += log(transIntronU);
    }
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
