#include "Sensor.NStart.h"

/*************************************************************
 **                       SensorNStart                      **
 *************************************************************/

extern Parameters PAR;

// ----------------------
// Default constructor.
// ----------------------
SensorNStart :: SensorNStart ()
{
  *Start = NULL;

  startP = PAR.getD("NStart.startP");
  startB = PAR.getD("NStart.startB");
}

// ----------------------
//  Default destructor.
// ----------------------
SensorNStart :: ~SensorNStart ()
{
  delete [] Start[0];
  delete [] Start[1];
}

// ----------------------
//  Init start.
// ----------------------
void SensorNStart :: Init (DNASeq *X)
{
  char tempname[FILENAME_MAX+1];

  type = Type_Start;

  if(*Start != NULL) {
    delete [] Start[0];
    delete [] Start[1];
  }

  Start[0] = new REAL[X->SeqLen+1];
  Start[1] = new REAL[X->SeqLen+1];
  
  fprintf(stderr, "Reading start file (NetStart).................");
  fflush(stderr);
  strcpy(tempname,PAR.getC("fstname"));
  strcat(tempname,".starts");
  ReadNStart(tempname, X->SeqLen, 0);
  fprintf(stderr,"forward,");
  fflush(stderr);
  
  strcpy(tempname,PAR.getC("fstname"));
  strcat(tempname,".startsR");
  ReadNStart(tempname, X->SeqLen, 1);
  fprintf(stderr," reverse done\n");
  
  CheckStart(X, Start);
}

// ----------------------
//  Read start file.
// ----------------------
void SensorNStart :: ReadNStart (char name[FILENAME_MAX+1], int Len, int rev)
{
  FILE *fp;
  int i,j = -1;
  double force;

  if (!(fp = fopen(name, "r"))) {
    fprintf(stderr, "cannot open start file %s\n", name);
    exit(2);
  }

  for (i = 0; i<= Len; i++)
    Start[rev][i] = 0.0;

  while (1) {
    i = fscanf(fp,"%d %lf %*s\n", &j, &force);
    if (i == EOF) break;
    if (i < 2) {
      fprintf(stderr, "Error in start file %s, position %d\n", name, j);
      exit(2);
    }
    if (rev) j = Len-j+2;
    Start[rev][j-1] = pow(force, startB)*(exp(-startP));
  }
  if (j == -1) fprintf(stderr,"WARNING: empty NetStart file !\n");
  fclose(fp);
}

// ------------------------
//  GiveInfo signal start.
// ------------------------
void SensorNStart :: GiveInfo (DNASeq *X, int pos, DATA *d)
{
  d->Start[0] += Start[0][pos];
  d->Start[1] += Start[1][pos];
}
