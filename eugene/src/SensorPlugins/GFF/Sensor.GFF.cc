#include "Sensor.GFF.h"

/*************************************************************
 **                        GFFObject                        **
 *************************************************************/
// -------------------------
//  Default constructor.
// -------------------------
GFFObject :: GFFObject ()
{
  feature[0] = '\000';
  start  = end   = 0;
  strand = frame = '?';
}

// -------------------------
//  Default constructor.
// -------------------------
GFFObject :: GFFObject (char* fea, int s, int e, char st, char fr)
{
  strcpy(feature,fea);
  start = s;
  end   = e;
  strand = st;
  frame  = fr;
}

// -------------------------
//  Print GFFObject
// -------------------------
void GFFObject :: Print ()
{
  printf("\t\t%s\t%d\t%d\t\t%c\t%c\n",feature, start, end, strand, frame);
}

// -------------------------
//  Print GFF Header
// -------------------------
void GFFObject :: PrintHeader ()
{
  printf("name\tsource\tfeature\tstart\tend\tscore\tstrand\tframe\n");
}

// -------------------------
//  Default destructor.
// -------------------------
GFFObject :: ~GFFObject () {}


/*************************************************************
 **                       SensorGFF                         **
 *************************************************************/

extern Parameters PAR;

// ----------------------
// Default constructor.
// ----------------------
SensorGFF :: SensorGFF (int n) : Sensor(n)
{
}

// ----------------------
//  Default destructor.
// ----------------------
SensorGFF :: ~SensorGFF ()
{
  // Clear the data structures
  gffList.clear();
}

// ----------------------
//  Init.
// ----------------------
void SensorGFF :: Init (DNASeq *X)
{
  char tempname[FILENAME_MAX+1];

  // Type initialisation
  type = Type_Multiple;
  
  // Clear the data structure
  gffList.clear();

  fprintf(stderr, "Reading GFF file.............");
  fflush(stderr);
  strcpy(tempname,PAR.getC("fstname"));
  strcat(tempname,".gff");
  ReadGFF(tempname);
  fprintf(stderr, "done\n");
  fflush(stderr);

  if (PAR.getI("Output.graph")) Plot(X);

  //   for (int i=0; i<(int)gffList.size(); i++) {
  //   if(i==0)
  //   gffList[i]->PrintHeader();
  //   gffList[i]->Print();
  //   fflush(stdout);
  //   }
}

// --------------------------
//  Read start forward file.
// --------------------------
void SensorGFF :: ReadGFF (char name[FILENAME_MAX+1])
{
  FILE *fp;
  int i;
  char *feature;
  int  start, end;
  char strand, frame;

  if (!(fp = fopen(name, "r"))) {
    fprintf(stderr, "cannot open GFF file %s\n", name);
    exit(2);
  }
  
  int j=0;
  while (1) {
    j++;
    i = fscanf(fp,"%*s %*s %s %d %d %*s %c %c",
	       feature, &start, &end, &strand, &frame);
    if (i < 5) {
      if (i==-1) {
	if(j==1)
	  fprintf(stderr,"WARNING: empty GFF file !...");
      }
      else {
	fprintf(stderr, "Error in GFF file %s, line %d.\n", name, j);
	exit(2);
      }
    }
    else {
      if (frame != '.')
	gffList.push_back(new GFFObject(feature, start, end, strand, frame));
    }
    char c = fgetc(fp);
    while (c != '\n' && c != EOF) c=fgetc(fp);
    if (fgetc(fp) == EOF) break;
  }
  fclose(fp);
}

// ------------------------
//  GiveInfo.
// ------------------------
void SensorGFF :: GiveInfo (DNASeq *X, int pos, DATA *d)
{
}

// ----------------------------
//  Plot Sensor information.
// ----------------------------
void SensorGFF :: Plot(DNASeq *X)
{
  const int PredWidth = 2;
  
  for(int i=0; i<(int)gffList.size(); i++) {
    int f = ((int)gffList[i]->frame - 48) + 1;
    if (gffList[i]->strand == '-')
      f *= -1;
          
    for (int j=gffList[i]->start; j<=gffList[i]->end; j++)
      PlotBarI(j, f, 0.27, PredWidth, 9);
    
    if (gffList[i]->strand == '.') {
      f *= -1;
      for (int j=gffList[i]->start; j<=gffList[i]->end; j++)
 	PlotBarI(j, f, 0.27, PredWidth, 9);
    }
  }
}

// ------------------
//  Post analyse.
// ------------------
void SensorGFF :: PostAnalyse(Prediction *pred)
{
}
