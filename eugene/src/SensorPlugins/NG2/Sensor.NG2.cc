#include "Sensor.NG2.h"

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
  
  iterAccF = iterAccR = iterDonF = iterDonR = 0;
  
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
  
  CheckSplices(X, vPosAccF, vPosDonF, vPosAccR, vPosDonR);
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

// -----------------------
//  ResetIter.
// -----------------------
void SensorNG2 :: ResetIter ()
{
  iterAccF = iterAccR = iterDonF = iterDonR = 0;
}

// ------------------------
//  GiveInfo signal NG2.
// ------------------------
void SensorNG2 :: GiveInfo (DNASeq *X, int pos, DATA *d)
{
  if( iterAccF <= (int)vPosAccF.size()  &&  vPosAccF[iterAccF] == pos ) {
    if(d->Acc[0] == 0.0)
      d->Acc[0] = vValAccF[iterAccF];
    iterAccF++;
  }
  if( iterAccR <= (int)vPosAccR.size()  &&  vPosAccR[iterAccR] == pos ) {
    if(d->Acc[1] == 0.0)
      d->Acc[1] = vValAccR[iterAccR];
    iterAccR++;
  }
  if( iterDonF <= (int)vPosDonF.size()  &&  vPosDonF[iterDonF] == pos ) {
    if(d->Don[0] == 0.0)
      d->Don[0] = vValDonF[iterDonF];
    iterDonF++;
  }
  if( iterDonR <= (int)vPosDonR.size()  &&  vPosDonR[iterDonR] == pos ) {
    if(d->Don[1] == 0.0)
      d->Don[1] = vValDonR[iterDonR];
    iterDonR++;
  }
}

// ------------------------
//  GiveInfoAt signal NG2.
// ------------------------
void SensorNG2 :: GiveInfoAt (DNASeq *X, int pos, DATA *d)
{
  iter = find(vPosAccF.begin(), vPosAccF.end(), pos);
  if(iter != vPosAccF.end() && d->Acc[0] == 0.0)
    d->Acc[0] = vValAccF[iter-vPosAccF.begin()];

  iter = find(vPosAccR.begin(), vPosAccR.end(), pos);
  if(iter != vPosAccR.end() && d->Acc[1] == 0.0)
    d->Acc[1] = vValAccR[iter-vPosAccR.begin()];

  iter = find(vPosDonF.begin(), vPosDonF.end(), pos);
  if(iter != vPosDonF.end() && d->Don[0] == 0.0)
    d->Don[0] = vValDonF[iter-vPosDonF.begin()];

  iter = find(vPosDonR.begin(), vPosDonR.end(), pos);
  if(iter != vPosDonR.end() && d->Don[0] == 0.0)
    d->Don[1] = vValDonR[iter-vPosDonR.begin()];
}
