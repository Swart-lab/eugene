#include "Sensor.Riken.h"

/*************************************************************
 **                        SensorRiken                      **
 *************************************************************/

extern Parameters PAR;

bool Before(const RAFLgene A, const RAFLgene B)
{ return (A.deb < B.deb); };

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
 
  RAFLgene tmp;

  std::vector <RAFLgene> RAFLtmp;
  
  fprintf(stderr, "Reading RAFL gene............");
  fflush(stderr);
  strcpy(tempname, PAR.getC("fstname"));
  strcat(tempname, ".riken");
  fRAFL = FileOpen(NULL, tempname, "r");
  while ((j=fscanf (fRAFL,"%d %d %*s %d %d %*s %s",
		    &beg5, &end5, &beg3, &end3, name)) == 5) {
    tmp.deb = Min(beg3, Min(end3, Min(beg5, end5)));
    tmp.fin = Max(beg3, Max(end3, Max(beg5, end5)));
    strcpy(tmp.ID, name);
    tmp.sens = (((beg5+end5) <= (beg3+end3)) ? 1 : -1);

    // si le gene est trop court, on ne peut pas connaitre le sens !
    if (abs((beg5+end5) - (beg3+end3)) <100) tmp.sens = 0;

    if (abs(tmp.deb-tmp.fin) > MAX_RIKEN_LENGTH) {
      fprintf(stderr, "\nWARNING: Check RAFL data: Riken %s rejected, too long transcript (%d kb)\n",
	      name, abs(tmp.deb-tmp.fin)/1000);
      continue;
    }
    if (abs(beg5-end5) > MAX_RIKEN_EST_LENGTH) {
      fprintf(stderr, "\nWARNING: Check RAFL data: Riken %s rejected, 5' EST mapped to %d bp\n",
	      name, abs(beg5-end5));
      continue;
    }
    if (abs(beg3-end3) > MAX_RIKEN_EST_LENGTH) {
      fprintf(stderr, "\nWARNING: Check RAFL data: Riken %s rejected, 3'' EST mapped to %d bp\n",
	      name, abs(beg3-end3));
      continue;
    }
    if (abs(tmp.deb-tmp.fin) < MIN_RIKEN_LENGTH) {
      fprintf(stderr, "\nWARNING: Check RAFL data: Riken %s rejected, too short transcript (%d bp)\n",
	      name, abs(tmp.deb-tmp.fin));
      continue;
      }
    if (abs(beg5-end5) < MIN_RIKEN_EST_LENGTH) {
      fprintf(stderr, "\nWARNING: Check RAFL data: Riken %s rejected, 5' EST mapped to %d bp\n",
	      name, abs(beg5-end5));
      continue;
    }
    if (abs(beg3-end3) < MIN_RIKEN_EST_LENGTH) {
      fprintf(stderr, "\nWARNING: Check RAFL data: Riken %s rejected, 3'' EST mapped to %d bp\n",
	      name, abs(beg3-end3));
      continue;
    }
    if (tmp.deb <= 0  || tmp.fin > X->SeqLen) {
      fprintf(stderr, "\nWARNING: Check RAFL data: Riken %s rejected, coordinate(s) %d-%d out of range\n",
      	      name, tmp.deb, tmp.fin);
      continue;
    }
    
    RAFLtmp.push_back(tmp);
  }

  if (j != EOF) {
    fprintf(stderr, "Incorrect RAFL file\n");
    exit(2);
  }

  fclose(fRAFL);
  fprintf(stderr, "%d RAFL EST pairs read, ", RAFLtmp.size());

  sort(RAFLtmp.begin(), RAFLtmp.end(), Before);

  int RAFLindice=0;
  int RAFLtmpindice=0;

  //  for (RAFLtmpindice=0; RAFLtmpindice< RAFLtmp.size(); RAFLtmpindice++) {
  //  fprintf(stderr, "\nRAFLtmp[%d]: deb=%d fin=%d sens=%d ID=%s\n",RAFLtmpindice,RAFLtmp[RAFLtmpindice].deb,RAFLtmp[RAFLtmpindice].fin,RAFLtmp[RAFLtmpindice].sens,RAFLtmp[RAFLtmpindice].ID);
  // }

    if ((int)RAFLtmp.size()>0) RAFL.push_back(RAFLtmp[RAFLtmpindice]); // on retient le premier

    for (RAFLtmpindice=1; RAFLtmpindice<(int)RAFLtmp.size(); RAFLtmpindice++) { 
      // on analyse chaque tmp que l'on compare avec le dernier rafl retenu
  
      if (RAFL[RAFLindice].fin < RAFLtmp[RAFLtmpindice].deb) { // pas d'overlap
	RAFL.push_back(RAFLtmp[RAFLtmpindice]);
	RAFLindice++;
	continue; // tt va bien, on passe au tmp suivant
      }
      if (RAFL[RAFLindice].fin - RAFLtmp[RAFLtmpindice].deb >= 60) { // cas de grand overlap
	fprintf(stderr,"fusion %s-%s...",RAFL[RAFLindice].ID,RAFLtmp[RAFLtmpindice].ID);
	if ( (RAFL[RAFLindice].sens!=0) && 
	     (RAFLtmp[RAFLtmpindice].sens!=0) && 
	     (RAFL[RAFLindice].sens!=RAFLtmp[RAFLtmpindice].sens) &&
	     (RAFL[RAFLindice].sens!=2) ) { // orientations contraires
	  fprintf(stderr, "\nWARNING: Check RAFL data: contig containing rikens in different orientations (%s and %s)\n",RAFL[RAFLindice].ID,RAFLtmp[RAFLtmpindice].ID);
	  RAFL[RAFLindice].sens = 2;
	}
	// grand overlap, fusion des deux riken (meme gene)
	RAFL[RAFLindice].deb = Min( RAFL[RAFLindice].deb , RAFLtmp[RAFLtmpindice].deb );
	RAFL[RAFLindice].fin = Max( RAFL[RAFLindice].fin , RAFLtmp[RAFLtmpindice].fin );
	if (RAFL[RAFLindice].sens==0) 
	  RAFL[RAFLindice].sens = RAFLtmp[RAFLtmpindice].sens;
	continue;
      }

      else { // little overlap, real overlapping genes
	fprintf(stderr,"%s-%s overlapping...",RAFL[RAFLindice].ID,RAFLtmp[RAFLtmpindice].ID);
	j= RAFLtmp[RAFLtmpindice].deb; // on memorise l'endroit ou tronquer le dernier rafl
	i=RAFLtmpindice; // on raccourcit le debut de tous les petits chevauchant:
	while ( (i<(int)RAFLtmp.size()) && (RAFL[RAFLindice].fin >= RAFLtmp[i].deb) ) {
	  RAFLtmp[i].deb= RAFL[RAFLindice].fin;
	  i++;
	}
	RAFL[RAFLindice].fin = j;// on raccourcit la fin du dernier
	
	RAFL.push_back(RAFLtmp[RAFLtmpindice]);
	RAFLindice++;
	continue;
      }
    }

    fprintf(stderr,"resulting %d\n",RAFL.size());
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
    if (pos > RAFL[RAFLindex].fin) {
      // si on depasse le RAFL, on prend l'eventuel prochain
      if (RAFLindex+1 < (int)RAFL.size()) RAFLindex++; 
      else RAFL_A_Traiter = FALSE;
    }

    // dans le Riken:
    if (pos >= RAFL[RAFLindex].deb-2) {
      // bordure min:
      if (pos == RAFL[RAFLindex].deb-2) {
	if ((RAFL[RAFLindex].sens== 1) || (RAFL[RAFLindex].sens== -1))
	  RAFLpos = ( (RAFL[RAFLindex].sens == 1) ? 1 : 2);
	else 
	  RAFLpos=0;
      }
      else {
	if (pos < RAFL[RAFLindex].fin) RAFLpos = 3;
	else {
	  // frontiere max:
	  if (pos == RAFL[RAFLindex].fin) {
	    if ((RAFL[RAFLindex].sens== 1) || (RAFL[RAFLindex].sens== -1))
	      RAFLpos = ( (RAFL[RAFLindex].sens == -1) ? 1 : 2);
	    else 
	      RAFLpos=0;
	  }
	}
      }
    }

    for(int i=0; i<3; i++)
      if ((RAFLpos == 1) || (RAFLpos == 2) ||
	  ((RAFLpos == 3) && (RAFL[RAFLindex].sens == -1)) )
	d->contents[i] += RAFLPenalty;   //ExonF

    for(int i=3; i<6; i++)
      if ((RAFLpos == 1) || (RAFLpos == 2) ||
	  ((RAFLpos == 3) && (RAFL[RAFLindex].sens == 1)) )
	d->contents[i] += RAFLPenalty;   //ExonR

    if ( (RAFLpos == 1) || (RAFLpos == 2) || 
	 ( (RAFLpos == 3) && (RAFL[RAFLindex].sens == -1)) )
      d->contents[6] += RAFLPenalty;     //IntronF

    if ( (RAFLpos == 1) || (RAFLpos == 2) ||
	 ( (RAFLpos == 3) && (RAFL[RAFLindex].sens == 1)) )
      d->contents[7] += RAFLPenalty;     //IntronR

    if (RAFLpos == 3) {
      d->contents[8] += RAFLPenalty;     //InterG
      
      // Warning SIGNAL START MODIFICATION //
      d->sig[DATA::Start].weight[Signal::ForwardNo] = 0.0;
      d->sig[DATA::Start].weight[Signal::ReverseNo] = 0.0;
    }
    
    if ( (RAFLpos == 2) || ((RAFL[RAFLindex].sens == -1) && (RAFLpos >= 1) ) )
      d->contents[9] += RAFLPenalty;     //UTR5'F
    
    if ( (RAFLpos == 2) || ((RAFL[RAFLindex].sens == 1) && (RAFLpos >= 1) ) )
      d->contents[10]+= RAFLPenalty;     //UTR5'R
    
    if ( (RAFLpos == 1) || ((RAFL[RAFLindex].sens == -1) && (RAFLpos >= 1) ) )
      d->contents[11]+= RAFLPenalty;     //UTR3'F
    
    if ( (RAFLpos == 1) || ((RAFL[RAFLindex].sens == 1) && (RAFLpos >= 1) ) )
      d->contents[12]+= RAFLPenalty;     //UTR3'R
  }
}

// ----------------------------
//  Plot Sensor information
// ----------------------------
void SensorRiken :: Plot(DNASeq *TheSeq)
{
  const int HLen = 30;
  int plotcolor=0;

  for (int j =0; j<(int)RAFL.size() ;j++) {
    plotcolor= ( (RAFL[j].sens==2) ? 2 : 2-RAFL[j].sens);
    PlotBarF(RAFL[j].deb, 0, 0.9, 0.2, plotcolor);
    PlotLine(RAFL[j].deb, RAFL[j].deb+HLen, 0, 0, 1.0, 1.0, plotcolor);
    PlotBarF(RAFL[j].fin, 0, 0.9, 0.2, plotcolor);
    PlotLine(RAFL[j].fin-HLen,RAFL[j].fin, 0, 0, 1.0, 1.0, plotcolor);
  }
}

// ------------------
//  Post analyse
// ------------------
void SensorRiken :: PostAnalyse(Prediction *pred)
{
}
