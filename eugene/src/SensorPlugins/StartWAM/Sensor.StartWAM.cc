#include "Sensor.StartWAM.h"
//#define NORM(x,n) (((n)+(Max(-(n),x)))/(n))
#define NORM(x,n) ( (Min(x,n))/(n))
#define plotscoreincrease 7

/*******************************************************
 **                  SensorStartWAM                   **
 *******************************************************/

extern Parameters PAR;

// ----------------------
//  Default constructor.
// ----------------------
SensorStartWAM :: SensorStartWAM (int n) : Sensor(n)
{
  int i;
  char* filename;
  type = Type_Start;
  char TPstartfile[FILENAME_MAX+1];
  char FPstartfile[FILENAME_MAX+1];

  strcpy(TPstartfile,PAR.getC("StartWAM.TPstartfile"));
  strcpy(FPstartfile,PAR.getC("StartWAM.FPstartfile"));
  order = PAR.getI("StartWAM.order");
  coef = PAR.getD("StartWAM.coef");
  pen = PAR.getD("StartWAM.pen");
  beforestart = PAR.getI("StartWAM.beforestart");
  afterstart = PAR.getI("StartWAM.afterstart");
  startsitelen = beforestart + 3 + afterstart;
  ADN = new ChaineADN;
  TPSTARTMOD = new TabChaine<ChaineADN,unsigned short int>[startsitelen](order,ADN);
  FPSTARTMOD = new TabChaine<ChaineADN,unsigned short int>[startsitelen](order,ADN);

// Loading models
// (=binary matrix files containing a markov model, one per position of each signal)
  fprintf (stderr,"Loading WAM models... Start ");
  fflush (stderr);
  for (i=0; i<startsitelen ;i++) {
    fprintf (stderr,"%d ",i);
    fflush (stderr);
    filename= new char[FILENAME_MAX+1];
    sprintf(filename,"%s",TPstartfile);
    if (i<10) sprintf(filename+strlen(PAR.getC("StartWAM.TPstartfile")),"0%d",i);
    else sprintf(filename+strlen(PAR.getC("StartWAM.TPstartfile")),"%d",i);
    ReadModelFile(TPSTARTMOD[i],filename);
    delete filename;
    filename= new char[FILENAME_MAX+1];
    sprintf(filename,"%s",FPstartfile);
    if (i<10) sprintf(filename+strlen(PAR.getC("StartWAM.FPstartfile")),"0%d",i);
    else sprintf(filename+strlen(PAR.getC("StartWAM.FPstartfile")),"%d",i);
    ReadModelFile(FPSTARTMOD[i],filename);
    delete filename;
  }

//  printf ("----------- MODELS AFFICHAGE ----------\n");
//  fflush (stdout);
//  for (int i=0; i<accsitelen ;i++) {
//    printf("Affichage TPACCMOD[%d]\n",i);
//    TPACCMOD[i].affichagevaleurs(2);
//  }
  fprintf (stderr,"... done\n");
  fflush (stderr);
}

// ----------------------
//  Default destructor.
// ----------------------
SensorStartWAM :: ~SensorStartWAM ()
{
  delete [] TPSTARTMOD;
  delete [] FPSTARTMOD;
  delete ADN;
}

// ---------------------
//  Init splice.
// ----------------------
void SensorStartWAM :: Init (DNASeq *X)
{
  if (PAR.getI("Output.graph")) Plot(X);
}

void SensorStartWAM :: 
ReadModelFile (TabChaine<ChaineADN, unsigned short int> &model, char* filename)
{
  FILE* fp;
  if (! (fp = FileOpen(EugDir, filename, "rb"))) {
    fprintf(stderr, "cannot open model file %s\n", filename);
    exit(2);
  }
  if (model.chargefichier(fp)) {
    fprintf(stderr,"Error when loading model file %s\n",filename);
    exit(2);
  }
  fclose (fp);
}

double SensorStartWAM :: 
startscoring (char* oligont,int order, int insitepos)
{
  int i;
  for (i=0; (unsigned)i<strlen(oligont);i++) {
    // test if an unknown nucleotide is present (like "n"), out of the alphabet
    if ((unsigned)ADN->operator[](oligont[i]) == ADN->taille)
      return 0.0; // don't talk
  }
  return (log(TPSTARTMOD[insitepos].usi2real(TPSTARTMOD[insitepos].proba(oligont,order))) -
	  log(FPSTARTMOD[insitepos].usi2real(FPSTARTMOD[insitepos].proba(oligont,order)))   );
  // likelihood ratio: log ( proba(nt with true model)/proba(nt with false model) )
}

// -----------------------
//  GiveInfo signal start.
// -----------------------
void SensorStartWAM :: GiveInfo (DNASeq *X, int pos, DATA *d)
{
  // i= seq X indice, j= oligont indice, insitepos = nt position in the site
  int i,j,insitepos;
  double score;
  char* oligont = new char[order+2];
  oligont[order+1] ='\0'; // oligo stores a short word for asking scoring functions

  ////////// START Forward (need enough context) //////////////
  if ( (X->IsEStart(pos,1)) &&
       (pos-beforestart-order > 0) && 
       (pos+2+afterstart < X->SeqLen) ) {
    score=0.0;
    insitepos=0;
    for (i= pos-beforestart-order; i<= pos+2+afterstart-order; i++) {
      for (j=0;j<=order;j++) {
	oligont[j] = toupper((*X)[i+j]);
      }
      score += startscoring (oligont,order,insitepos);
      insitepos++;
    }
    d->sig[DATA::Start].weight[Signal::Forward] += (score * coef) - pen;
//    d->sig[DATA::Start].weight[Signal::ForwardNo] -= score;
  }

  ////////// START Reverse (need enough context)   //////////////

  if ( (X->IsEStart(pos-1,-1)) &&
       (pos-1+beforestart+order < X->SeqLen) &&
       (pos-3-afterstart > 0) ) {
    score=0.0;
    insitepos=0;
    for (i= pos-1+beforestart+order; i >= pos-3-afterstart+order; i--) {
      for (j=0;j<=order;j++) {
	oligont[j] = toupper((*X)(i-j));
      }
      score += startscoring (oligont,order,insitepos);
      insitepos++;
    }
    d->sig[DATA::Start].weight[Signal::Reverse] += (score * coef) - pen;
//    d->sig[DATA::Start].weight[Signal::ReverseNo] -= score;
  }
}
// ----------------------------
//  Plot Sensor information
// ----------------------------
void SensorStartWAM :: Plot(DNASeq *X)
{  
  int pos;
  DATA data;

  for(pos=0; pos < X-> SeqLen ; pos++){

    data.sig[DATA::Start].weight[Signal::Forward] = 0.0 ;
    data.sig[DATA::Start].weight[Signal::Reverse] = 0.0 ;

    GiveInfo (X, pos, &data);

    if (data.sig[DATA::Start].weight[Signal::Forward] != 0) {
      data.sig[DATA::Start].weight[Signal::Forward] += plotscoreincrease;
      if (data.sig[DATA::Start].weight[Signal::Forward] > 0)
	PlotBarF(pos,(pos%3)+1,0.5,NORM(data.sig[DATA::Start].weight[Signal::Forward],10.0),2);
    }
    if (data.sig[DATA::Start].weight[Signal::Reverse] != 0) {
      data.sig[DATA::Start].weight[Signal::Reverse] += plotscoreincrease;
      if (data.sig[DATA::Start].weight[Signal::Reverse] > 0)
	PlotBarF(pos,-((X->SeqLen-pos-1)%3)-1,0.5,NORM(data.sig[DATA::Start].weight[Signal::Reverse],10.0),2);
    }
  }
}

// ------------------
//  Post analyse
// ------------------
void SensorStartWAM :: PostAnalyse(Prediction *pred)
{
}
