#include "Sensor.Riken.h"

/*************************************************************
 **                        SensorRiken                      **
 *************************************************************/

extern Parameters PAR;

bool Before(const RAFLgene *A, const RAFLgene *B)
{ return (A->deb < B->deb); };

// ************
// * RAFLgene *     // RAFL: Riken Arabidopsis Full Length cDNA 
// ************
// ----------------------
// Default constructor.
// ----------------------
RAFLgene :: RAFLgene () { deb = fin = sens = 0; ID[0] = '0'; }

// ----------------------
// Default destructor.
// ----------------------
RAFLgene :: ~RAFLgene () {}

// ****************
// * Sensor Riken *
// ****************
// ----------------------
// Default constructor.
// ----------------------
SensorRiken :: SensorRiken (int n) : Sensor(n)
{
}

// ----------------------
//  Default destructor.
// ----------------------
SensorRiken :: ~SensorRiken ()
{
  RAFL.clear();
}

// ----------------------
//  Init Riken.
// ----------------------
void SensorRiken :: Init (DNASeq *X)
{
  RAFLpos = RAFLindex = 0;
  RAFL_A_Traiter = TRUE;
  FILE *fRAFL;
  int i, j;
  int beg5, end5;
  int beg3, end3;
  char name[FILENAME_MAX+1];
  char tempname[FILENAME_MAX+1];

  type = Type_Content;
 
  RAFL.clear();

  RAFLgene *tmp = new RAFLgene();
  std::vector <RAFLgene*> RAFLtmp;
  
  fprintf(stderr, "Reading RAFL gene............");
  fflush(stderr);
  strcpy(tempname, PAR.getC("fstname"));
  strcat(tempname, ".riken");
  fRAFL = FileOpen(NULL, tempname, "r");
  while ((j=fscanf (fRAFL,"%d %d %*s %d %d %*s %s",
		    &beg5, &end5, &beg3, &end3, name)) == 5) {
    tmp->deb = Min(beg3, Min(end3, Min(beg5, end5)));
    tmp->fin = Max(beg3, Max(end3, Max(beg5, end5)));
    strcpy(tmp->ID, name);
    
    tmp->sens = (((beg5+end5) < (beg3+end3)) ? 1 : -1);
    // si le gene est trop court, on ne peut pas connaitre le sens !
    if (abs((beg5+end5) - (beg3+end3)) <100) tmp->sens = 0;
    RAFLtmp.push_back(tmp);
  }
  fclose(fRAFL);
  if (j != EOF) {
    fprintf(stderr, "Incorrect RAFL file\n");
    exit(2);
  }
  fprintf(stderr, "%d RAFL EST pairs read, ", RAFLtmp.size());

  sort(RAFLtmp.begin(), RAFLtmp.end(), Before);
  
  for (j=0; j<(int)RAFLtmp.size()-1; j++) {
    if (RAFLtmp[j]->fin - RAFLtmp[j+1]->deb >= 60) { // grand overlap
      fprintf(stderr, "fusion...");
      // grand overlap, on fusionne les deux riken (meme gene)
      RAFLtmp[j]->deb = Min( RAFLtmp[j]->deb, RAFLtmp[j+1]->deb );
      RAFLtmp[j]->fin = Max( RAFLtmp[j]->fin, RAFLtmp[j+1]->fin );
      if ( (RAFLtmp[j]->sens == RAFLtmp[j+1]->sens) || 
	   ((RAFLtmp[j]->sens==0) || (RAFLtmp[j+1]->sens==0))) {
	RAFLtmp[j]->sens=((RAFLtmp[j]->sens==0) ? 
			  RAFLtmp[j+1]->sens : RAFLtmp[j]->sens);
	RAFL.push_back(RAFLtmp[j]); 
	// si pas de contradiction dans les sens, on prend le gene resultant 
      }
      j++;
      //on saute le suivant
    }
    else{
      if (RAFLtmp[j]->fin - RAFLtmp[j+1]->deb > 0) {
	fprintf(stderr, "overlap...");
	i = RAFLtmp[j]->fin;
	RAFLtmp[j]->fin   = Min(i, RAFLtmp[j+1]->deb);
	RAFLtmp[j+1]->deb = Max(i, RAFLtmp[j+1]->deb);
      }
      RAFL.push_back(RAFLtmp[j]);
    }
  }
  if(j==(int)RAFLtmp.size()-1) RAFL.push_back(RAFLtmp[j]);
  
  fprintf(stderr,"%d kept\n",RAFL.size());
  fflush(stderr);
  if (RAFL.size() < 1) RAFL_A_Traiter = FALSE;

  if (PAR.getI("Output.graph")) Plot(X);
}

// -------------------------
//  GiveInfo Content Riken.
// -------------------------
void SensorRiken :: GiveInfo (DNASeq *X, int pos, DATA *d)
{
  // Calcul de la position par rapport aux genes RAFL (Riken Ara.Full-Length)
  // valeurs de RAFLpos:
  // 0-> en dehors,
  // 1-> frontiere 5prime (intergenique ou UTR5  obligatoire),
  // 2-> frontiere 3prime (intergenique ou UTR3  obligatoire),
  // 3-> dedans(penalisation IG)

  if(RAFL_A_Traiter) {
    RAFLpos = 0;
    if (pos > RAFL[RAFLindex]->fin) {
      // si on depasse le RAFL, on prend l'eventuel prochain
      if (RAFLindex+1 < (int)RAFL.size()) RAFLindex++; 
      else RAFL_A_Traiter = FALSE;
    }
    
    // dans le Riken:
    if (pos >= RAFL[RAFLindex]->deb-2) {
      // bordure min:
      if (pos == RAFL[RAFLindex]->deb-2)
	RAFLpos = ( (RAFL[RAFLindex]->sens == 1) ? 1 : 2);
      else {
	if (pos < RAFL[RAFLindex]->fin) RAFLpos = 3;
	else {
	  // frontiere max:
	  if (pos == RAFL[RAFLindex]->fin) {
	    RAFLpos = ( (RAFL[RAFLindex]->sens == -1) ? 1 : 2);
	  }
	}
      }
    }
    
    for(int i=0; i<3; i++)
      if ((RAFLpos == 1) || (RAFLpos == 2) ||
	  ((RAFLpos == 3) && (RAFL[RAFLindex]->sens == -1)) )
	d->ContentScore[i] += RAFLPenalty;   //ExonF
    
    for(int i=3; i<6; i++)
      if ((RAFLpos == 1) || (RAFLpos == 2) ||
	  ((RAFLpos == 3) && (RAFL[RAFLindex]->sens == 1)) )
	d->ContentScore[i] += RAFLPenalty;   //ExonR
    
    if ( (RAFLpos == 1) || (RAFLpos == 2) || 
	 ( (RAFLpos == 3) && (RAFL[RAFLindex]->sens == -1)) )
      d->ContentScore[6] += RAFLPenalty;     //IntronF
    
    if ( (RAFLpos == 1) || (RAFLpos == 2) ||
	 ( (RAFLpos == 3) && (RAFL[RAFLindex]->sens == 1)) )
      d->ContentScore[7] += RAFLPenalty;     //IntronR
    
    if (RAFLpos == 3)
      d->ContentScore[8] += RAFLPenalty;     //InterG
    
    if ( (RAFLpos == 2) || ((RAFL[RAFLindex]->sens == -1) && (RAFLpos >= 1) ) )
      d->ContentScore[9] += RAFLPenalty;     //UTR5'F
    
    if ( (RAFLpos == 2) || ((RAFL[RAFLindex]->sens == 1) && (RAFLpos >= 1) ) )
      d->ContentScore[10]+= RAFLPenalty;     //UTR5'R
    
    if ( (RAFLpos == 1) || ((RAFL[RAFLindex]->sens == -1) && (RAFLpos >= 1) ) )
      d->ContentScore[11]+= RAFLPenalty;     //UTR3'F
    
    if ( (RAFLpos == 1) || ((RAFL[RAFLindex]->sens == 1) && (RAFLpos >= 1) ) )
      d->ContentScore[12]+= RAFLPenalty;     //UTR3'R
  }
}

// ----------------------------
//  Plot Sensor information
// ----------------------------
void SensorRiken :: Plot(DNASeq *TheSeq)
{
  const int HLen = 30;
  
  for (int j =0; j<(int)RAFL.size() ;j++) {
    PlotBarF(RAFL[j]->deb, 0, 0.9, 0.2, 2-RAFL[j]->sens);
    PlotLine(RAFL[j]->deb, RAFL[j]->deb+HLen, 0, 0, 1.0, 1.0, 2-RAFL[j]->sens);
    PlotBarF(RAFL[j]->fin, 0, 0.9, 0.2, 2-RAFL[j]->sens);
    PlotLine(RAFL[j]->fin-HLen,RAFL[j]->fin, 0, 0, 1.0, 1.0, 2-RAFL[j]->sens);
  }
}

// ------------------
//  Post analyse
// ------------------
void SensorRiken :: PostAnalyse(Prediction *pred)
{
}
