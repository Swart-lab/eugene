#include "Sensor.GSplicer.h"

/*************************************************************
 **                       SensorGSplicer
 *************************************************************/

extern Parameters PAR;

// ----------------------
// Default constructor.
// ----------------------
SensorGSplicer :: SensorGSplicer (int n) : Sensor(n)
{
  coefAcc = PAR.getD("GSplicer.coefAcc");
  penAcc  = PAR.getD("GSplicer.penAcc");
  coefDon = PAR.getD("GSplicer.coefDon");
  penDon  = PAR.getD("GSplicer.penDon");
}

// ----------------------
//  Default destructor.
// ----------------------
SensorGSplicer :: ~SensorGSplicer ()
{
  // Clear the data structures 
  vPosAccF.clear();  vPosAccR.clear();
  vPosDonF.clear();  vPosDonR.clear();
  vValAccF.clear();  vValAccR.clear();
  vValDonF.clear();  vValDonR.clear();
}

// ----------------------
//  Init.
// ----------------------
void SensorGSplicer :: Init (DNASeq *X)
{
  char tempname[FILENAME_MAX+1];

  // Type initialisation
  type = Type_Splice;

  iAccF = iDonF = 0;
  
  // Clear the data structures 
  vPosAccF.clear();  vPosAccR.clear();
  vPosDonF.clear();  vPosDonR.clear();
  vValAccF.clear();  vValAccR.clear();
  vValDonF.clear();  vValDonR.clear();

  fprintf(stderr, "Reading splice site file (GeneSplicer)........");
  fflush(stderr);
  strcpy(tempname,PAR.getC("fstname"));
  strcat(tempname,".Gsplicer");
  fprintf(stderr,"forward, reverse ");
  fflush(stderr);
  ReadGSplicer(tempname, X->SeqLen);
  fprintf(stderr,"done\n");
  
  CheckSplices(X,vPosAccF, vPosDonF, vPosAccR, vPosDonR);
  
  if (PAR.getI("Output.graph")) Plot(X);
  
  iAccR = (int)vPosAccR.size() - 1;
  iDonR = (int)vPosDonR.size() - 1;
}

// -------------------------------------
//  Read GeneSplicer file.
// -------------------------------------
void SensorGSplicer :: ReadGSplicer(char name[FILENAME_MAX+1], int SeqLen)
{
  FILE *fp;
  char type[10];
  int i,j,nt1,nt2;
  double score;
  
  if (!(fp = fopen(name, "r"))) {
    fprintf(stderr, "cannot open splice sites file %s\n",  name);
    exit(2);
  }
  
  i=0;
  while(!feof(fp)){
    i++;
    j = fscanf(fp,"%d %d %lf %*s %s\n", &nt1, &nt2, &score, type);
    
    if (j < 4) {
      fprintf(stderr, "\nError in splice sites file! %s, line %d\n", name, i);
      exit(2);
    }
        
    // Forward
    if ( (nt2-nt1) == 1 ) {
      if(type[0]=='d') {
	vPosDonF.push_back( nt1-1 );
	vValDonF.push_back( score );
      }
      else {
	vPosAccF.push_back( nt2 );
	vValAccF.push_back( score );
      }
    }
    // Reverse
    if ( (nt2-nt1) == -1 ) {
      if(type[0]=='d') {
	vPosDonR.push_back( nt1 );
	vValDonR.push_back( score );
      }
      else {
	vPosAccR.push_back( nt2-1 );
	vValAccR.push_back( score );
      }
    }
  }
  fclose(fp);
}

// ------------------------
//  GiveInfo.
// ------------------------
void SensorGSplicer :: GiveInfo (DNASeq *X, int pos, DATA *d)
{
  // Accepteur Forward
  if((iAccF != 0                    &&  vPosAccF[iAccF-1] >= pos) ||
     (iAccF < (int)vPosAccF.size()  &&  vPosAccF[iAccF]   <  pos))
    {
      iter = lower_bound(vPosAccF.begin(), vPosAccF.end(), pos);
      if(*iter == pos) {
	d->sig[DATA::Acc].weight[Signal::Forward] +=
	  (vValAccF[iter-vPosAccF.begin()] * coefAcc) - penAcc;
	//d->sig[DATA::Acc].weight[Signal::ForwardNo] += ...
	iAccF = iter-vPosAccF.begin() + 1;
      }
      else iAccF = iter-vPosAccF.begin();
    }
  else if(iAccF < (int)vPosAccF.size()  &&  vPosAccF[iAccF] == pos)
    {
      d->sig[DATA::Acc].weight[Signal::Forward] += (vValAccF[iAccF] * coefAcc) - penAcc;
      //d->sig[DATA::Acc].weight[Signal::ForwardNo] += ...
      iAccF++;
    }
  
  // Accepteur Reverse
  if((iAccR < (int)vPosAccR.size()  &&  vPosAccR[iAccR+1] >= pos) ||
     (iAccR > -1                    &&  vPosAccR[iAccR]   <  pos))
    {
      iter = lower_bound(vPosAccR.begin(), vPosAccR.end(), pos, std::greater<int>());
      if(*iter == pos) {
	d->sig[DATA::Acc].weight[Signal::Reverse] +=
	  (vValAccR[iter-vPosAccR.begin()] * coefAcc) - penAcc;
	//d->sig[DATA::Acc].weight[Signal::ReverseNo] += ...
	iAccR = iter-vPosAccR.begin();
      }
      else iAccR = iter-vPosAccR.begin() - 1;
    }
  else if(iAccR > -1  &&  vPosAccR[iAccR] == pos)
    {
      d->sig[DATA::Acc].weight[Signal::Reverse] += (vValAccR[iAccR] * coefAcc) - penAcc;
      //d->sig[DATA::Acc].weight[Signal::ReverseNo] += ...
      iAccR--;
    }
  
  // Donneur Forward
  if((iDonF != 0                    &&  vPosDonF[iDonF-1] >= pos) ||
     (iDonF < (int)vPosDonF.size()  &&  vPosDonF[iDonF]   <  pos))
    {
      iter = lower_bound(vPosDonF.begin(), vPosDonF.end(), pos);
      if(*iter == pos) {
	d->sig[DATA::Don].weight[Signal::Forward] +=
	  (vValDonF[iter-vPosDonF.begin()] * coefDon) - penDon;
	//d->sig[DATA::Don].weight[Signal::ForwardNo] += ...
	iDonF = iter-vPosDonF.begin() + 1;
      }
      else iDonF = iter-vPosDonF.begin();
    }
  else if(iDonF < (int)vPosDonF.size()  &&  vPosDonF[iDonF] == pos)
    {
      d->sig[DATA::Don].weight[Signal::Forward] += (vValDonF[iDonF] * coefDon) - penDon;
      //d->sig[DATA::Don].weight[Signal::ForwardNo] += ...
      iDonF++;
    }
  
  // Donneur Reverse
  if((iDonR < (int)vPosDonR.size()  &&  vPosDonR[iDonR+1] >= pos) ||
     (iDonR > -1                    &&  vPosDonR[iDonR]   <  pos))
    {
      iter = lower_bound(vPosDonR.begin(), vPosDonR.end(), pos, std::greater<int>());
      if(*iter == pos) {
	d->sig[DATA::Don].weight[Signal::Reverse] += 
	  (vValDonR[iter-vPosDonR.begin()] * coefDon) - penDon;
	//d->sig[DATA::Don].weight[Signal::ReverseNo] += ...
	iDonR = iter-vPosDonR.begin();
      }
      else iDonR = iter-vPosDonR.begin() - 1;
    }
  else if(iDonR > -1  &&  vPosDonR[iDonR] == pos)
    {
      d->sig[DATA::Don].weight[Signal::Reverse] += (vValDonR[iDonR] * coefDon) - penDon;
      //d->sig[DATA::Don].weight[Signal::ReverseNo] += ...
      iDonR--;
    }
}

// ----------------------------
//  Plot Sensor information.
// ----------------------------
void SensorGSplicer :: Plot(DNASeq *X)
{
}

// ------------------
//  Post analyse.
// ------------------
void SensorGSplicer :: PostAnalyse(Prediction *pred)
{
}
