#include "Sensor.NStart.h"
#define NORM(x,n) (((n)+(Max(-(n),x)))/(n))

/*************************************************************
 **                       SensorNStart                      **
 *************************************************************/

extern Parameters PAR;

// ----------------------
// Default constructor.
// ----------------------
SensorNStart :: SensorNStart (int n) : Sensor(n)
{
  startP = PAR.getD("NStart.startP");
  startB = PAR.getD("NStart.startB");

  PositionGiveInfo = -1;
}

// ----------------------
//  Default destructor.
// ----------------------
SensorNStart :: ~SensorNStart ()
{
  vPosF.clear();
  vValF.clear();
  vPosR.clear();
  vValR.clear();
}

// ----------------------
//  Init start.
// ----------------------
void SensorNStart :: Init (DNASeq *X)
{
  char tempname[FILENAME_MAX+1];

  type = Type_Start;
  
  indexR = 0;
  
  vPosF.clear();
  vValF.clear();
  vPosR.clear();
  vValR.clear();
  
  fprintf(stderr, "Reading start file (NetStart).................");
  fflush(stderr);
  strcpy(tempname,PAR.getC("fstname"));
  strcat(tempname,".starts");
  ReadNStartF(tempname, X->SeqLen);
  fprintf(stderr,"forward,");
  fflush(stderr);
  
  strcpy(tempname,PAR.getC("fstname"));
  strcat(tempname,".startsR");
  ReadNStartR(tempname, X->SeqLen);
  fprintf(stderr," reverse done\n");
  
  CheckStart(X,vPosF, vPosR);

  if (PAR.getI("Output.graph")) Plot(X);

  indexF = (int)vPosR.size() - 1;
}

// --------------------------
//  Read start forward file.
// --------------------------
void SensorNStart :: ReadNStartF (char name[FILENAME_MAX+1], int Len)
{
  FILE *fp;
  int i,j = -1;
  double force;

  if (!(fp = fopen(name, "r"))) {
    fprintf(stderr, "cannot open start file %s\n", name);
    exit(2);
  }
  
  while (1) {
    i = fscanf(fp,"%d %lf %*s\n", &j, &force);
    if (i == EOF) break;
    if (i < 2) {
      fprintf(stderr, "Error in start file %s, position %d\n", name, j);
      exit(2);
    }
    vPosF.push_back( j-1 );
    vValF.push_back( pow(force, startB)*(exp(-startP)) );
  }
  if (j == -1) fprintf(stderr,"WARNING: empty NetStart file !\n");
  fclose(fp);
}

// --------------------------
//  Read start reverse file.
// --------------------------
void SensorNStart :: ReadNStartR (char name[FILENAME_MAX+1], int Len)
{
  FILE *fp;
  int i,j = -1;
  double force;

  if (!(fp = fopen(name, "r"))) {
    fprintf(stderr, "cannot open start file %s\n", name);
    exit(2);
  }

  while (1) {
    i = fscanf(fp,"%d %lf %*s\n", &j, &force);
    if (i == EOF) break;
    if (i < 2) {
      fprintf(stderr, "Error in start file %s, position %d\n", name, j);
      exit(2);
    }
    j = Len-j+2;
    vPosR.push_back( j-1 );
    vValR.push_back( pow(force, startB)*(exp(-startP)) );
  }
  if (j == -1) fprintf(stderr,"WARNING: empty NetStart file !\n");
  fclose(fp);
}

// ------------------------
//  GiveInfo signal start.
// ------------------------
void SensorNStart :: GiveInfo (DNASeq *X, int pos, DATA *d)
{
  bool update = false;

  if ( (PositionGiveInfo == -1) || (pos != PositionGiveInfo+1) ) update = true; // update indexes on vectors
  PositionGiveInfo = pos;
  
  // Start Forward
  if(!vPosF.empty()) {
    if (update) 
      indexF = lower_bound(vPosF.begin(), vPosF.end(), pos)-vPosF.begin();
    
    if((indexF<(int)vPosF.size()) && (vPosF[indexF] == pos)) {
      d->sig[DATA::Start].weight[Signal::Forward] += log(vValF[indexF]);
      d->sig[DATA::Start].weight[Signal::ForwardNo] += log(1.0-vValF[indexF]);
      indexF++;
    }
  }
  
  // Start Reverse
  if (!vPosR.empty()) {
    if (update) { 
      indexR = lower_bound(vPosR.begin(), vPosR.end(), pos, std::greater<int>())-vPosR.begin();
      if (indexR==(int)vPosR.size()) indexR--; // if pos is before first site, then point first site
      if (vPosR[indexR]<pos) indexR = -1; // if pos is after last site, then do not point
    }

    if((indexR!=-1) && (indexR<(int)vPosR.size()) && (vPosR[indexR] == pos)) {
      d->sig[DATA::Start].weight[Signal::Reverse] += log(vValR[indexR]);
      d->sig[DATA::Start].weight[Signal::ReverseNo] += log(1.0-vValR[indexR]);
      indexR--;
    }
  }

}

// ----------------------------
//  Plot Sensor information
// ----------------------------
void SensorNStart :: Plot(DNASeq *X)
{
  for (int i =0; i < (int)vPosF.size(); i++)
    PlotBarF(vPosF[i],(vPosF[i]%3)+1,0.5,NORM(log(vValF[i]),4.0),2);

  for (int i =0; i < (int)vPosR.size(); i++)
    PlotBarF(vPosR[i],-((X->SeqLen-vPosR[i])%3)-1,0.5,NORM(log(vValR[i]),4.0),2);
}

// ------------------
//  Post analyse
// ------------------
void SensorNStart :: PostAnalyse(Prediction *pred)
{
}
