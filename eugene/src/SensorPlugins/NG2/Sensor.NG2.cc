#include "Sensor.NG2.h"
#define NORM(x,n) (((n)+(Max(-(n),x)))/(n))

/*************************************************************
 **                        SensorNetGene2                   **
 *************************************************************/

extern Parameters PAR;

// ----------------------
//  Default constructor.
// ----------------------
SensorNG2 :: SensorNG2 (int n) : Sensor(n)
{
  accB = PAR.getD("NG2.accB",GetNumber());
  accP = PAR.getD("NG2.accP",GetNumber());
  donB = PAR.getD("NG2.donB",GetNumber());
  donP = PAR.getD("NG2.donP",GetNumber());

  PositionGiveInfo = -1;
}

// ----------------------
//  Default destructor.
// ----------------------
SensorNG2 :: ~SensorNG2 ()
{
  vPosAccF.clear();  vPosAccR.clear();
  vPosDonF.clear();  vPosDonR.clear();
  vValAccF.clear();  vValAccR.clear();
  vValDonF.clear();  vValDonR.clear();
}

// ----------------------
//  Init NG2.
// ----------------------
void SensorNG2 :: Init (DNASeq *X)
{
  char tempname[FILENAME_MAX+1];

  type = Type_Splice;
  
  iAccF = iDonF = 0;
  
  vPosAccF.clear();  vPosAccR.clear();
  vPosDonF.clear();  vPosDonR.clear();
  vValAccF.clear();  vValAccR.clear();
  vValDonF.clear();  vValDonR.clear();
  
  fprintf(stderr, "Reading splice site file (NetGene2)...........");  
  fflush(stderr);
  strcpy(tempname,PAR.getC("fstname"));
  strcat(tempname,".splices");
  ReadNG2F(tempname, X->SeqLen);
  fprintf(stderr,"forward,");
  fflush(stderr);

  strcpy(tempname,PAR.getC("fstname"));
  strcat(tempname,".splicesR");
  ReadNG2R(tempname, X->SeqLen);
  fprintf(stderr," reverse done\n");
  
  CheckSplices(X,vPosAccF, vPosDonF, vPosAccR, vPosDonR);
  
  if (PAR.getI("Output.graph")) Plot(X);
  
  iAccR = (int)vPosAccR.size() - 1;
  iDonR = (int)vPosDonR.size() - 1;
}

// -----------------------------
//  Read NetGene2 forward file.
// -----------------------------
void SensorNG2 :: ReadNG2F(char name[FILENAME_MAX+1], int SeqLen)
{
  FILE *fp;
  char buf[FILENAME_MAX];
  char sacc[10],sdon[10];
  char altsacc[10],altsdon[10];
  int i = 0, j;
  
  if (!(fp = fopen(name, "r"))) {
    fprintf(stderr, "cannot open splice sites file %s\n", name);
    exit(2);
  }
  fgets(buf,FILENAME_MAX-1,fp);
  
  for (i = 0; i < SeqLen; i++) {
    // final value of netgene2 
    j = fscanf(fp,"%*s %*s %s %s %*s %*s %*s %*s %s %s %*s %*s",
	       altsdon,altsacc,sdon,sacc);
    if (j < 4) {
      fprintf(stderr, "Error in splice sites file %s, line %d\n", name, i+2);
      exit(2);
     }
    if (sdon[0] == '-') strcpy(sdon,altsdon);
    if (sacc[0] == '-') strcpy(sacc,altsacc);
    
    if( atof(sacc) != 0.0 ) {
      vPosAccF.push_back( i+1 );
      vValAccF.push_back( pow(atof(sacc),  accB)*accP);
    }
    if( atof(sdon) != 0.0 ) {
      vPosDonF.push_back( i );
      vValDonF.push_back( pow(atof(sdon),  donB)*donP);
    }
  }
  fclose(fp);
}

// -----------------------------
//  Read NetGene2 reverse file.
// -----------------------------
void SensorNG2 :: ReadNG2R(char name[FILENAME_MAX+1], int SeqLen)
{
  FILE *fp;
  char buf[FILENAME_MAX];
  char sacc[10],sdon[10];
  char altsacc[10],altsdon[10];
  int i = 0, j, k;
  
  if (!(fp = fopen(name, "r"))) {
    fprintf(stderr, "cannot open splice sites file %s\n", name);
    exit(2);
  }
  fgets(buf,FILENAME_MAX-1,fp);

  k = SeqLen;
  for (i = 0; i < SeqLen; i++) {
    // final value of netgene2 
    j = fscanf(fp,"%*s %*s %s %s %*s %*s %*s %*s %s %s %*s %*s",
	       altsdon,altsacc,sdon,sacc);
    if (j < 4) {
      fprintf(stderr, "Error in splice sites file %s, line %d\n", name, i+2);
      exit(2);
     }
    if (sdon[0] == '-') strcpy(sdon,altsdon);
    if (sacc[0] == '-') strcpy(sacc,altsacc);

    if( atof(sacc) != 0.0 ) {
      vPosAccR.push_back( k-1 );
      vValAccR.push_back( pow(atof(sacc),  accB)*accP);
    }
    if( atof(sdon) != 0.0 ) {
      vPosDonR.push_back( k );
      vValDonR.push_back( pow(atof(sdon),  donB)*donP);
    }
    k--;
  }
  fclose(fp);
}

// ------------------------
//  GiveInfo signal NG2.
// ------------------------
void SensorNG2 :: GiveInfo (DNASeq *X, int pos, DATA *d)
{
  bool update = false;

  if ( (PositionGiveInfo == -1) || (pos != PositionGiveInfo+1) ) update = true; // update indexes on vectors
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
void SensorNG2 :: Plot(DNASeq *X)
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
void SensorNG2 :: PostAnalyse(Prediction *pred)
{
}
