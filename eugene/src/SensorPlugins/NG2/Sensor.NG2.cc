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
  accB = PAR.getD("NG2.accB");
  accP = PAR.getD("NG2.accP");
  donB = PAR.getD("NG2.donB");
  donP = PAR.getD("NG2.donP");
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
