#include "SensorNG2.h"

/*************************************************************
 **                        SensorNetGene2                   **
 *************************************************************/

extern Parameters PAR;

// ----------------------
//  Default constructor.
// ----------------------
SensorNG2 :: SensorNG2 ()
{
  *Acc = NULL;
}

// ----------------------
//  Default destructor.
// ----------------------
SensorNG2 :: ~SensorNG2 ()
{
  delete [] Acc[0];
  delete [] Acc[1];
  delete [] Don[0];
  delete [] Don[1];
}

// ----------------------
//  Init NG2.
// ----------------------
void SensorNG2 :: Init (DNASeq *X)
{
  char tempname[FILENAME_MAX+1];

  type = Type_Splice;

  if(*Acc != NULL) {
    delete [] Acc[0];
    delete [] Acc[1];
    delete [] Don[0];
    delete [] Don[1];
  }

  Acc[0] = new REAL[X->SeqLen+1];
  Don[0] = new REAL[X->SeqLen+1];
  Acc[1] = new REAL[X->SeqLen+1];
  Don[1] = new REAL[X->SeqLen+1];
  for(int i=0; i<=X->SeqLen; i++)
    Acc[0][i] = Don[0][i] = Acc[1][i] = Don[1][i] = 0.0; 

  fprintf(stderr, "Reading splice site file (NetGene2)...........");  
  fflush(stderr);
  strcpy(tempname,PAR.fstname);
  strcat(tempname,".splices");
  ReadNG2(tempname, 1, X->SeqLen);
  fprintf(stderr,"forward,");
  fflush(stderr);

  strcpy(tempname,PAR.fstname);
  strcat(tempname,".splicesR");
  ReadNG2(tempname, -1, X->SeqLen);
  fprintf(stderr," reverse done\n");
  
  CheckSplices(X, Acc, Don);
}

// ----------------------
// Read NetGene2 file.
// ----------------------
void SensorNG2 :: ReadNG2(char name[FILENAME_MAX+1], int dir, int SeqLen)
{
  FILE *fp;
  char buf[FILENAME_MAX];
  char sacc[10],sdon[10];
  char altsacc[10],altsdon[10];
  int i = 0, j, k;
  int rev;

  if (dir == 1) rev = 0; else rev = 1;

  if (!(fp = fopen(name, "r"))) {
    fprintf(stderr, "cannot open splice sites file %s\n", name);
    exit(2);
  }
  fgets(buf,FILENAME_MAX-1,fp);

  if (dir == 1) k = 0; else k = SeqLen;
  Acc[rev][i] = 0.0;
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

    Acc[rev][k+dir] = atof(sacc);
    Acc[rev][k+dir] = pow(Acc[rev][k+dir],PAR.AccB[0])*PAR.AccP[0];
    
    Don[rev][k] = atof(sdon);
    Don[rev][k] = pow(Don[rev][k],PAR.DonB[0])*PAR.DonP[0];
    
    k += dir;
  }
  Don[rev][k] = 0.0;  
  fclose(fp);
}

// ------------------------
//  GiveInfo signal NG2.
// ------------------------
void SensorNG2 :: GiveInfo (DNASeq *X, int pos, DATA *d)
{
  if(d->Acc[0] == 0.0) d->Acc[0] = Acc[0][pos];
  if(d->Acc[1] == 0.0) d->Acc[1] = Acc[1][pos];
  if(d->Don[0] == 0.0) d->Don[0] = Don[0][pos];
  if(d->Don[1] == 0.0) d->Don[1] = Don[1][pos];
}
