#include "Sensor.Repeat.h"

/*************************************************************
 **                        SensorRepeat                     **
 *************************************************************/

extern Parameters PAR;

// ----------------------
// Default constructor.
// ----------------------
SensorRepeat :: SensorRepeat (int n) : Sensor(n)
{
}

// ----------------------
//  Default destructor.
// ----------------------
SensorRepeat :: ~SensorRepeat ()
{
  vPos.clear();
}

// ----------------------
//  Init Repeat.
// ----------------------
void SensorRepeat :: Init (DNASeq *X)
{
  char tempname[FILENAME_MAX+1];
  FILE* ncfile;
  int i, deb, fin;

  type = Type_Content;

  index = 0;

  vPos.clear();
  
  fprintf(stderr,"Reading Intergenic regions... ");
  fflush(stderr);
  
  strcpy(tempname,PAR.getC("fstname"));
  strcat(tempname,".ig");
  ncfile = FileOpen(NULL,tempname, "r");
  
  while (fscanf(ncfile,"%d %d\n", &deb, &fin) != EOF)  {
    deb = Max(1,deb)-1;
    fin = Min(X->SeqLen,fin)-1;
    for (i = deb; i <= fin; i++) {
      vPos.push_back( i );
      if (PAR.getI("Output.graph")) PlotBarI(i,0,0.25,2,6);
    }
  }
}

// -----------------------
//  ResetIter.
// -----------------------
void SensorRepeat :: ResetIter ()
{
  index = 0;
}

// --------------------------
//  GiveInfo Content Repeat.
// --------------------------
void SensorRepeat :: GiveInfo (DNASeq *X, int pos, DATA *d)
{
  if( index <= (int)vPos.size()  &&  vPos[index] == pos ) {
    for(int i=0; i<8; i++)   // Exon(6) + Intron(2)
      d->ContentScore[i] += IGPenalty;
    index++;
  }
}

// ----------------------------
//  GiveInfoAt Content Repeat.
// ----------------------------
void SensorRepeat :: GiveInfoAt (DNASeq *X, int pos, DATA *d)
{
  iter = find(vPos.begin(), vPos.end(), pos);
  if(iter != vPos.end())
    for(int i=0; i<8; i++)   // Exon(6) + Intron(2)
      d->ContentScore[i] += IGPenalty;
}
