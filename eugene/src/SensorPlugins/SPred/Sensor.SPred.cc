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
  PositionGiveInfo = -1;
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

  // vectors for reverse are put in the increasing order
  reverse(vPosAccR.begin(), vPosAccR.end()); 
  reverse(vValAccR.begin(), vValAccR.end());
  reverse(vPosDonR.begin(), vPosDonR.end());
  reverse(vValDonR.begin(), vValDonR.end());

  iAccF = iDonF = iAccR = iDonR = 0;

  if (PAR.getI("Output.graph")) Plot(X);
}

// -------------------------------------
//  Read Splice Predictor forward file.
// -------------------------------------
void SensorSPred :: ReadSPredF(char name[FILENAME_MAX+1], int SeqLen)
{
  FILE *fp;
  char buf[FILENAME_MAX];
  double par[3],strength;
  char *type;
  int i=1,j, k;
  int prevkA = -1, prevkD = -1;
  
  fp = FileOpen(NULL,name,"r");
  
  k = -1;
  while (1) {

    if (!fgets(buf,FILENAME_MAX-1,fp))
      break;
    
    // vieux format ou -p5 ?
    if (buf[0] == 'A' || buf[0] == 'D') { //old
      j = sscanf(buf+9,"%d",&k);
      j += sscanf(buf+38,"%lf",&par[0]);
      j += sscanf(buf+45,"%lf",&par[1]);
      j += sscanf(buf+52,"%lf",&par[2]);
      type = buf;
    }
    else
      {
	j = sscanf(buf,"%*d %*s %*s %d %*f %*f %*f %lf %lf %lf",&k, par,par+1,par+2);
	type = buf+2;
      }

    strength = par[0];

    // erreur: on a pas tout lu ou ca ne croit pas
      if ((j < 4) || ((type[0] == 'D') && (k < prevkD)) || ((k < prevkA)))      {
	fprintf(stderr,"%d %c %d %d %d\n",j,type[0],k,prevkD,prevkA);
 	fprintf(stderr, "\nError in splice sites file %s, line %d\n", name, i);
 	exit(2);
      }
      
      if (type[0] == 'D' && strength != 0.0 && k>prevkD) {
	prevkD = k;
	vPosDonF.push_back( k-1 );
	vValDonF.push_back( pow(strength, donB)*donP );
      }  else
	if (strength != 0.0 && k > prevkA)  {
	  prevkA = k;
	  vPosAccF.push_back( k );
	  vValAccF.push_back( pow(strength, accB)*accP );
	}
      i++;
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
  double strength;
  double par[3];
  char* type;
  int i = 1, j, k;
  int prevkA = INT_MAX, prevkD = INT_MAX;

  fp = FileOpen(NULL,name,"r");

  k = -1;
  while (1) {
    if (!fgets(buf,FILENAME_MAX-1,fp))
      break;

    // vieux format ou -p5 ?
    if (buf[0] == 'A' || buf[0] == 'D') { //old
      j = sscanf(buf+9,"%d",&k);
      j += sscanf(buf+38,"%lf",&par[0]);
      j += sscanf(buf+45,"%lf",&par[1]);
      j += sscanf(buf+52,"%lf",&par[2]);
      type = buf;
    }
    else
      {
	j = sscanf(buf,"%*d %*s %*s %d %*f %*f %*f %lf %lf %lf",&k, par,par+1,par+2);
	type = buf+2;
      }

    strength = par[0];

    // on ne lit pas tout on ca l'index position ne decroit pas
    if ((j < 4) || (type[0] == 'D' && k > prevkD) || (k > prevkA)) {
      fprintf(stderr, "\nError in splice sites file %s, line %d\n", name, i);
      exit(2);
    }
          
    if (strength != 0.0) 
      if (type[0] == 'D') {
	prevkD = k;
	vPosDonR.push_back( k );
	vValDonR.push_back( pow(strength, donB)*donP );
      }
      else {
	prevkA = k;
	vPosAccR.push_back( k-1 );
	vValAccR.push_back( pow(strength, accB)*accP);
      }
    i++;
  }
  
  if (k == -1) fprintf(stderr,"WARNING: Empty splice predictor file !\n");
  fclose(fp);
}

// ------------------------
//  GiveInfo signal SPred.
// ------------------------
void SensorSPred :: GiveInfo (DNASeq *X,int pos, DATA *d)
{
  bool update = false;

  if ( (PositionGiveInfo == -1) || (pos != PositionGiveInfo+1) ) update = true; 
  PositionGiveInfo = pos;

  // Accepteur Forward
  if(!vPosAccF.empty()) {
    if (update) 
      iAccF = lower_bound(vPosAccF.begin(), vPosAccF.end(), pos)-vPosAccF.begin();
    
    if((iAccF<(int)vPosAccF.size()) && (vPosAccF[iAccF] == pos)) {
      d->sig[DATA::Acc].weight[Signal::Forward] += log(vValAccF[iAccF]);
      d->sig[DATA::Acc].weight[Signal::ForwardNo] += log(1.0-vValAccF[iAccF]);
      iAccF++;
    }
  }
  
  // Accepteur Reverse
  if (!vPosAccR.empty()) {
    if (update) 
      iAccR = lower_bound(vPosAccR.begin(), vPosAccR.end(), pos)-vPosAccR.begin();

    if((iAccR<(int)vPosAccR.size()) && (vPosAccR[iAccR] == pos)) {
      d->sig[DATA::Acc].weight[Signal::Reverse] += log(vValAccR[iAccR]);
      d->sig[DATA::Acc].weight[Signal::ReverseNo] += log(1.0-vValAccR[iAccR]);
      iAccR++;
    }
  }

  // Donneur Forward
  if(!vPosDonF.empty()) {
    if (update) 
      iDonF = lower_bound(vPosDonF.begin(), vPosDonF.end(), pos)-vPosDonF.begin();

    if ((iDonF<(int)vPosDonF.size()) && (vPosDonF[iDonF] == pos)) {
      d->sig[DATA::Don].weight[Signal::Forward] += log(vValDonF[iDonF]);
      d->sig[DATA::Don].weight[Signal::ForwardNo] += log(1.0-vValDonF[iDonF]);
      iDonF++;
    }
  }
  
  // Donneur Reverse
  if(!vPosDonR.empty()) {
    if (update) 
      iDonR = lower_bound(vPosDonR.begin(), vPosDonR.end(), pos)-vPosDonR.begin();

    if((iDonR<(int)vPosDonR.size()) && (vPosDonR[iDonR] == pos)) {
      d->sig[DATA::Don].weight[Signal::Reverse] += log(vValDonR[iDonR]);
      d->sig[DATA::Don].weight[Signal::ReverseNo] += log(1.0-vValDonR[iDonR]);
      iDonR++;
    }
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
