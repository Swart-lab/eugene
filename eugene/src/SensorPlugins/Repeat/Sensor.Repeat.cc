#include "SensorRepeat.h"

/*************************************************************
 **                        SensorRepeat                     **
 *************************************************************/

extern Parameters PAR;

// ----------------------
// Default constructor.
// ----------------------
SensorRepeat :: SensorRepeat ()
{
  ForcedIG = NULL;
}

// ----------------------
//  Default destructor.
// ----------------------
SensorRepeat :: ~SensorRepeat ()
{
  delete [] ForcedIG;
}

// ----------------------
//  Init Repeat.
// ----------------------
void SensorRepeat :: Init (DNASeq *X)
{
  char tempname[FILENAME_MAX+1];
  FILE* ncfile;
  int i, j, deb, fin;

  type = Type_Content;

  if(ForcedIG != NULL)
    delete [] ForcedIG;
  ForcedIG = NULL;
  
  fprintf(stderr,"Reading Intergenic regions... ");
  fflush(stderr);
  
  strcpy(tempname,PAR.getC("fstname"));
  strcat(tempname,".ig");
  ncfile = FileOpen(NULL,tempname, "r");
  
  ForcedIG = new unsigned char[X->SeqLen+1];
  for (j = 0; j <= X->SeqLen; j++) ForcedIG[j] = FALSE;
  
  while (fscanf(ncfile,"%d %d\n", &deb, &fin) != EOF)  {
    deb = Max(1,deb)-1;
    fin = Min(X->SeqLen,fin)-1;
    for (i = deb; i <= fin; i++) {
       ForcedIG[i] = TRUE;
       if (PAR.getI("graph")) PlotBarI(i,0,0.25,2,6);
    }
  }
}

// --------------------------
//  GiveInfo Content Repeat.
// --------------------------
void SensorRepeat :: GiveInfo(DNASeq *X, int pos, DATA *d)
{
  for(int i=0; i<8; i++)   // Exon(6) + Intron(2)
    if(ForcedIG[pos]) d->ContentScore[i] += IGPenalty;
}
