#include "Sensor.SPred.h"

/*************************************************************
 **                     SensorSplicePredictor               **
 *************************************************************/

extern Parameters PAR;

// ----------------------
//  Default constructor.
// ----------------------
SensorSPred :: SensorSPred ()
{
  *Acc = NULL;

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
  delete [] Acc[0];
  delete [] Acc[1];
  delete [] Don[0];
  delete [] Don[1];
}

// ----------------------
//  Init SPred.
// ----------------------
void SensorSPred :: Init (DNASeq *X)
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

  fprintf(stderr, "Reading splice site file (Splice Predictor)...");  
  fflush(stderr);
  strcpy(tempname,PAR.getC("fstname"));
  strcat(tempname,".spliceP");
  ReadSPred(tempname, 1, X->SeqLen);
  fprintf(stderr,"forward,");
  fflush(stderr);

  strcpy(tempname,PAR.getC("fstname"));
  strcat(tempname,".splicePR");
  ReadSPred(tempname, -1, X->SeqLen);
  fprintf(stderr," reverse done\n");

  CheckSplices(X, Acc, Don);
}

// -----------------------------
//  Read Splice Predictor file.
// ----------------------------- 
void SensorSPred :: ReadSPred(char name[FILENAME_MAX+1], int dir, int SeqLen)
{
  FILE *fp;
  char buf[FILENAME_MAX];
  double strength,add;
  int i = 0, j, k;
  int rev;
  
  if (dir == 1) rev = 0; else rev = 1;

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
        
    if (buf[0] == 'D') {
      if (Don[rev][k-(dir+1)/2] == 0.0)
	Don[rev][k-(dir+1)/2] = pow(strength, donB)*donP;
    }
    else
      if (Acc[rev][k+(dir-1)/2] == 0.0)
	Acc[rev][k+(dir-1)/2] = pow(strength, accB)*accP;
  }
  
  if (k == -1) fprintf(stderr,"WARNING: Empty splice predictor file !\n");
  fclose(fp);
}

// ------------------------
//  GiveInfo signal SPred.
// ------------------------
void SensorSPred :: GiveInfo (DNASeq *X, int pos, DATA *d)
{
  if(d->Acc[0] == 0.0) d->Acc[0] = Acc[0][pos];
  if(d->Acc[1] == 0.0) d->Acc[1] = Acc[1][pos];
  if(d->Don[0] == 0.0) d->Don[0] = Don[0][pos];
  if(d->Don[1] == 0.0) d->Don[1] = Don[1][pos];
}
