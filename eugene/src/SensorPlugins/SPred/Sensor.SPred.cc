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
  int i = 1, j, k;
  int prevkA = -1, prevkD = -1;
  
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
        
    // erreur: on a pas tout lu ou ca ne croit pas
    if ((j < 4) || ((buf[0] == 'D') && (k <= prevkD)) || ((k <= prevkA)))
      {
 	fprintf(stderr, "\nError in splice sites file %s, line %d\n", name, i);
 	exit(2);
      }
  
    if (buf[0] == 'D' && strength != 0.0) {
      prevkD = k;
      vPosDonF.push_back( k-1 );
      vValDonF.push_back( pow(strength, donB)*donP );
    }
    else
      if (strength != 0.0) {
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
  double strength,add;
  int i = 1, j, k;
  int prevkA = INT_MAX, prevkD = INT_MAX;

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
        
    // on ne lit pas tout on ca l'index position ne decroit pas
    if ((j < 4) || (buf[0] == 'D' && k >= prevkD) || (k >= prevkA)) {
      fprintf(stderr, "\nError in splice sites file %s, line %d\n", name, i);
      exit(2);
    }
          
    if (buf[0] == 'D' && strength != 0.0) {
      prevkD = k;
      vPosDonR.push_back( k );
      vValDonR.push_back( pow(strength, donB)*donP );
    }
    else
      if (strength != 0.0) {
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
  // update indexes on vectors
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
    if (update) { 
      iAccR = lower_bound(vPosAccR.begin(), vPosAccR.end(), pos, std::greater<int>())-vPosAccR.begin();
      if (iAccR==(int)vPosAccR.size()) iAccR--; // if pos is before first site, then point first site
      if (vPosAccR[iAccR]<pos) iAccR = -1; // if pos is after last site, then do not point
    }

    if((iAccR!=-1) && (iAccR<(int)vPosAccR.size()) && (vPosAccR[iAccR] == pos)) {
      d->sig[DATA::Acc].weight[Signal::Reverse] += log(vValAccR[iAccR]);
      d->sig[DATA::Acc].weight[Signal::ReverseNo] += log(1.0-vValAccR[iAccR]);
      iAccR--;
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
    if (update) {
      iDonR = lower_bound(vPosDonR.begin(), vPosDonR.end(), pos, std::greater<int>())-vPosDonR.begin();
      if (iDonR==(int)vPosDonR.size()) iDonR--;  // if pos is before first site, then point first site
      if (vPosDonR[iDonR]<pos) iDonR = -1; // if pos is after last site, then do not point
    }

    if((iDonR!=-1) && (vPosDonR[iDonR] == pos)) {
      d->sig[DATA::Don].weight[Signal::Reverse] += log(vValDonR[iDonR]);
      d->sig[DATA::Don].weight[Signal::ReverseNo] += log(1.0-vValDonR[iDonR]);
      iDonR--;
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
