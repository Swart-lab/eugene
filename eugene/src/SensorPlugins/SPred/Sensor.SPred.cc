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
  accP = PAR.getD("SPred.accP");
  accB = PAR.getD("SPred.accB");
  donP = PAR.getD("SPred.donP");
  donB = PAR.getD("SPred.donB");
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

  iterAccF = iterAccR = iterDonF = iterDonR = 0;
  
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

  CheckSplices(X, vPosAccF, vPosDonF, vPosAccR, vPosDonR);

  if (PAR.getI("Output.graph")) Plot(X);
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

// -----------------------
//  ResetIter.
// -----------------------
void SensorSPred :: ResetIter ()
{
  iterAccF = iterAccR = iterDonF = iterDonR = 0;
}

// ------------------------
//  GiveInfo signal SPred.
// ------------------------
void SensorSPred :: GiveInfo (DNASeq *X, int pos, DATA *d)
{
  int i;
  if( iterAccF < (int)vPosAccF.size()  &&  vPosAccF[iterAccF] == pos ) {
    d->Acc[0] = vValAccF[iterAccF];
    iterAccF++;
  }
  i = (int)vPosAccR.size();
  if( abs(iterAccR) < i  &&  vPosAccR[iterAccR + i-1] == pos ) {
    d->Acc[1] = vValAccR[iterAccR + i-1];
    iterAccR--;
  }
  if( iterDonF < (int)vPosDonF.size()  &&  vPosDonF[iterDonF] == pos ) {
    d->Don[0] = vValDonF[iterDonF];
    iterDonF++;
  }
  i = (int)vPosDonR.size();
  if( abs(iterDonR) < i  &&  vPosDonR[iterDonR + i-1] == pos ) {
    d->Don[1] = vValDonR[iterDonR + i-1];
    iterDonR--;
  }
}

// --------------------------
//  GiveInfoAt signal SPred.
// --------------------------
void SensorSPred :: GiveInfoAt (DNASeq *X, int pos, DATA *d)
{
  iter = lower_bound(vPosAccF.begin(), vPosAccF.end(), pos);
  if(*iter == pos)
    d->Acc[0] = vValAccF[iter-vPosAccF.begin()];
  
  iter = lower_bound(vPosAccR.begin(), vPosAccR.end(), pos, greater<int>());
  if(*iter == pos)
    d->Acc[1] = vValAccR[iter-vPosAccR.begin()];
  
  iter = lower_bound(vPosDonF.begin(), vPosDonF.end(), pos);
  if(*iter == pos)
    d->Don[0] = vValDonF[iter-vPosDonF.begin()];
  
  iter = lower_bound(vPosDonR.begin(), vPosDonR.end(), pos, greater<int>());
  if(*iter == pos)
    d->Don[1] = vValDonR[iter-vPosDonR.begin()];
}

// ----------------------------
//  Plot Sensor information
// ----------------------------
void SensorSPred :: Plot(DNASeq *X)
{
  for (int i =0; i < (int)vPosAccF.size(); i++)
    PlotBarF(vPosAccF[i],4,0.5,NORM(log(vValAccF[i]),20.0),4);
  
  for (int i =0; i < (int)vPosDonF.size(); i++)
    PlotBarF(vPosDonF[i],4,0.5,NORM(log(vValDonF[i]),20.0),5);
  
  for (int i =0; i < (int)vPosAccR.size(); i++)
    PlotBarF(vPosAccR[i],-4,0.5, NORM(log(vValAccR[i]),20.0),4);

  for (int i =0; i < (int)vPosDonR.size(); i++)
    PlotBarF(vPosDonR[i],-4,0.5,NORM(log(vValDonR[i]),20.0),5);
}

// ------------------
//  Post analyse
// ------------------
void SensorSPred :: PostAnalyse(Prediction *pred)
{
}
