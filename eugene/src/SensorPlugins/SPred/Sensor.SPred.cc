#include "Sensor.SPred.h"
#define NORM(x,n) (((n)+(Max(-(n),x)))/(n))

/*************************************************************
 **                     SensorSplicePredictor               **
 *************************************************************/

extern Parameters PAR;

// ----------------------
//  Default constructor.
// ----------------------
SensorSPred :: SensorSPred (int n) : Sensor(n)
{
  accP = PAR.getD("SPred.accP",GetNumber());
  accB = PAR.getD("SPred.accB",GetNumber());
  donP = PAR.getD("SPred.donP",GetNumber());
  donB = PAR.getD("SPred.donB",GetNumber());
}

// ----------------------
//  Default destructor.
// ----------------------
SensorSPred :: ~SensorSPred ()
{
  vPosAccF.clear();  vPosAccR.clear();
  vPosDonF.clear();  vPosDonR.clear();
  vValAccF.clear();  vValAccR.clear();
  vValDonF.clear();  vValDonR.clear();
}

// ----------------------
//  Init SPred.
// ----------------------
void SensorSPred :: Init (DNASeq *X)
{
  char tempname[FILENAME_MAX+1];

  type = Type_Splice;

  iAccF = iDonF = 0;
  
  vPosAccF.clear();  vPosAccR.clear();
  vPosDonF.clear();  vPosDonR.clear();
  vValAccF.clear();  vValAccR.clear();
  vValDonF.clear();  vValDonR.clear();
  
  fprintf(stderr, "Reading splice site file (Splice Predictor)...");  
  fflush(stderr);
  strcpy(tempname,PAR.getC("fstname"));
  strcat(tempname,".spliceP");
  ReadSPredF(tempname, X->SeqLen);
  fprintf(stderr,"forward,");
  fflush(stderr);

  strcpy(tempname,PAR.getC("fstname"));
  strcat(tempname,".splicePR");
  ReadSPredR(tempname, X->SeqLen);
  fprintf(stderr," reverse done\n");

  CheckSplices(X,vPosAccF, vPosDonF, vPosAccR, vPosDonR);

  if (PAR.getI("Output.graph")) Plot(X);
  
  iAccR = (int)vPosAccR.size() - 1;
  iDonR = (int)vPosDonR.size() - 1;
}

// -------------------------------------
//  Read Splice Predictor forward file.
// -------------------------------------
void SensorSPred :: ReadSPredF(char name[FILENAME_MAX+1], int SeqLen)
{
  FILE *fp;
  char buf[FILENAME_MAX];
  double strength,add;
  int i = 0, j, k;
  
  if (!(fp = fopen(name, "r"))) {
    fprintf(stderr, "cannot open splice sites file %s\n", name);
    exit(2);
  }
  
  k = -1;
  while (1) {
    if (!fgets(buf,99,fp))
      break;
    
    j = 0;
    j += sscanf(buf+9,"%d",&k);
    j += sscanf(buf+38,"%lf",&strength);
    j += sscanf(buf+45,"%lf",&add);
    j += sscanf(buf+52,"%lf",&add);
        
    if (j < 4) {
      fprintf(stderr, "Error in splice sites file %s, line %d\n", name, i+2);
      exit(2);
    }
        
    if (buf[0] == 'D' && strength != 0.0) {
      vPosDonF.push_back( k-1 );
      vValDonF.push_back( pow(strength, donB)*donP );
    }
    else
      if (strength != 0.0) {
	vPosAccF.push_back( k );
	vValAccF.push_back( pow(strength, accB)*accP );
      }
  }
  
  if (k == -1) fprintf(stderr,"WARNING: Empty splice predictor file !\n");
  fclose(fp);
}

// -------------------------------------
//  Read Splice Predictor reverse file.
// -------------------------------------
void SensorSPred :: ReadSPredR(char name[FILENAME_MAX+1], int SeqLen)
{
  FILE *fp;
  char buf[FILENAME_MAX];
  double strength,add;
  int i = 0, j, k;
  
  if (!(fp = fopen(name, "r"))) {
    fprintf(stderr, "cannot open splice sites file %s\n", name);
    exit(2);
  }

  k = -1;
  while (1) {
    if (!fgets(buf,99,fp))
      break;
    
    j = 0;
    j += sscanf(buf+9,"%d",&k);
    j += sscanf(buf+38,"%lf",&strength);
    j += sscanf(buf+45,"%lf",&add);
    j += sscanf(buf+52,"%lf",&add);
        
    if (j < 4) {
      fprintf(stderr, "Error in splice sites file %s, line %d\n", name, i+2);
      exit(2);
    }
          
    if (buf[0] == 'D' && strength != 0.0) {
      vPosDonR.push_back( k );
      vValDonR.push_back( pow(strength, donB)*donP );
    }
    else
      if (strength != 0.0) {
	vPosAccR.push_back( k-1 );
	vValAccR.push_back( pow(strength, accB)*accP);
      }
  }
  
  if (k == -1) fprintf(stderr,"WARNING: Empty splice predictor file !\n");
  fclose(fp);
}

// ------------------------
//  GiveInfo signal SPred.
// ------------------------
void SensorSPred :: GiveInfo (DNASeq *X,int pos, DATA *d)
{
  // Accepteur Forward
  if((iAccF != 0                    &&  vPosAccF[iAccF-1] >= pos) ||
     (iAccF < (int)vPosAccF.size()  &&  vPosAccF[iAccF]   <  pos))
    {
      iter = lower_bound(vPosAccF.begin(), vPosAccF.end(), pos);
      if(*iter == pos) {
	d->sig[DATA::Acc].weight[Signal::Forward] += log(vValAccF[iter-vPosAccF.begin()]);
	d->sig[DATA::Acc].weight[Signal::ForwardNo] += log(1.0-vValAccF[iter-vPosAccF.begin()]);
	iAccF = iter-vPosAccF.begin() + 1;
      }
      else iAccF = iter-vPosAccF.begin();
    }
  else if(iAccF < (int)vPosAccF.size()  &&  vPosAccF[iAccF] == pos)
    {
      d->sig[DATA::Acc].weight[Signal::Forward] += log(vValAccF[iAccF]);
      d->sig[DATA::Acc].weight[Signal::ForwardNo] += log(1.0-vValAccF[iAccF]);
      iAccF++;
    }
  
  // Accepteur Reverse
  if((iAccR < (int)vPosAccR.size()  &&  vPosAccR[iAccR+1] >= pos) ||
     (iAccR > -1                    &&  vPosAccR[iAccR]   <  pos))
    {
      iter = lower_bound(vPosAccR.begin(), vPosAccR.end(), pos, std::greater<int>());
      if(*iter == pos) {
	d->sig[DATA::Acc].weight[Signal::Reverse] += log(vValAccR[iter-vPosAccR.begin()]);
	d->sig[DATA::Acc].weight[Signal::ReverseNo] += log(1.0-vValAccR[iter-vPosAccR.begin()]);
	iAccR = iter-vPosAccR.begin();
      }
      else iAccR = iter-vPosAccR.begin() - 1;
    }
  else if(iAccR > -1  &&  vPosAccR[iAccR] == pos)
    {
      d->sig[DATA::Acc].weight[Signal::Reverse] += log(vValAccR[iAccR]);
      d->sig[DATA::Acc].weight[Signal::ReverseNo] += log(1.0-vValAccR[iAccR]);
      iAccR--;
    }
  
  // Donneur Forward
  if((iDonF != 0                    &&  vPosDonF[iDonF-1] >= pos) ||
     (iDonF < (int)vPosDonF.size()  &&  vPosDonF[iDonF]   <  pos))
    {
      iter = lower_bound(vPosDonF.begin(), vPosDonF.end(), pos);
      if(*iter == pos) {
	d->sig[DATA::Don].weight[Signal::Forward] += log(vValDonF[iter-vPosDonF.begin()]);
	d->sig[DATA::Don].weight[Signal::ForwardNo] += log(1.0-vValDonF[iter-vPosDonF.begin()]);
	iDonF = iter-vPosDonF.begin() + 1;
      }
      else iDonF = iter-vPosDonF.begin();
    }
  else if(iDonF < (int)vPosDonF.size()  &&  vPosDonF[iDonF] == pos)
    {
      d->sig[DATA::Don].weight[Signal::Forward] += log(vValDonF[iDonF]);
      d->sig[DATA::Don].weight[Signal::ForwardNo] += log(1.0-vValDonF[iDonF]);
      iDonF++;
    }
  
  // Donneur Reverse
  if((iDonR < (int)vPosDonR.size()  &&  vPosDonR[iDonR+1] >= pos) ||
     (iDonR > -1                    &&  vPosDonR[iDonR]   <  pos))
    {
      iter = lower_bound(vPosDonR.begin(), vPosDonR.end(), pos, std::greater<int>());
      if(*iter == pos) {
	d->sig[DATA::Don].weight[Signal::Reverse] += log(vValDonR[iter-vPosDonR.begin()]);
	d->sig[DATA::Don].weight[Signal::ReverseNo] += log(1.0-vValDonR[iter-vPosDonR.begin()]);
	iDonR = iter-vPosDonR.begin();
      }
      else iDonR = iter-vPosDonR.begin() - 1;
    }
  else if(iDonR > -1  &&  vPosDonR[iDonR] == pos)
    {
      d->sig[DATA::Don].weight[Signal::Reverse] += log(vValDonR[iDonR]);
      d->sig[DATA::Don].weight[Signal::ReverseNo] += log(1.0-vValDonR[iDonR]);
      iDonR--;
    }
}

// ----------------------------
//  Plot Sensor information
// ----------------------------
void SensorSPred :: Plot(DNASeq *X)
{
  for (int i =0; i < (int)vPosAccF.size(); i++)
    PlotBarF(vPosAccF[i],4,0.5,NORM(log(vValAccF[i]),20.0),4);
  
  for (int i =0; i < (int)vPosDonF.size(); i++)
    PlotBarF(vPosDonF[i],4,0.5,NORM(log(vValDonF[i]),20.0),11);
  
  for (int i =0; i < (int)vPosAccR.size(); i++)
    PlotBarF(vPosAccR[i],-4,0.5, NORM(log(vValAccR[i]),20.0),4);

  for (int i =0; i < (int)vPosDonR.size(); i++)
    PlotBarF(vPosDonR[i],-4,0.5,NORM(log(vValDonR[i]),20.0),11);
}

// ------------------
//  Post analyse
// ------------------
void SensorSPred :: PostAnalyse(Prediction *pred)
{
}
