#include "Sensor.BlastX.h"

/******************************************************************
 **                          SensorBlastX                        **
 ******************************************************************
 * Blastx against protein databases 1st simple version, we simply *
 * enhance coding probability according to phase Dangerous (NR    *
 * contains translation of frameshifted sequences) another        *
 * possibility would be too forbid introns/intergenic states      *
 * 8 levels of confidence may be used.                            *
 ******************************************************************/

extern Parameters PAR;

// ----------------------
// Default constructor.
// ----------------------
SensorBlastX :: SensorBlastX ()
{
  ProtMatch = NULL;

  keyBXLevel[0] = PAR.getD("BlastX.level0");
  keyBXLevel[1] = PAR.getD("BlastX.level1");
  keyBXLevel[2] = PAR.getD("BlastX.level2");
  keyBXLevel[3] = PAR.getD("BlastX.level3");
  keyBXLevel[4] = PAR.getD("BlastX.level4");
  keyBXLevel[5] = PAR.getD("BlastX.level5");
  keyBXLevel[6] = PAR.getD("BlastX.level6");
  keyBXLevel[7] = PAR.getD("BlastX.level7");
  minL8 = PAR.getI("EuGene.minL8");
}

// ----------------------
//  Default destructor.
// ----------------------
SensorBlastX :: ~SensorBlastX ()
{
  delete [] ProtMatch;
  delete [] ProtMatchLevel;
  delete [] ProtMatchPhase;
}

// ----------------------
//  Init blastX.
// ----------------------
void SensorBlastX :: Init (DNASeq *X)
{
  int i, j, k;
  int Len = X->SeqLen;
  const int LevelColor[3] = {6,7,8}; 
  FILE *fblast;
  int   overlap, deb, fin, phase, ProtDeb, ProtFin, PProtFin, level;
  float score;
  int   Pfin = 0, PPhase = 0;
  char  A[128], B[128];
  char *ProtId, *PProtId, *tmp;
  REAL GlobalScore;
  REAL PGlobalScore = 0;
  const int MaxOverlap = 10; 
  const int MaxHitLen  = 15000;
  char tempname[FILENAME_MAX+1];

  type = Type_Content;

  if(ProtMatch != NULL) {
    delete [] ProtMatch;
    delete [] ProtMatchLevel;
    delete [] ProtMatchPhase;
  }

  ProtMatch      = new REAL[Len+1];
  ProtMatchLevel = new REAL[Len+1];
  ProtMatchPhase = new int[Len+1];
  for (i = 0; i<= Len; i++) {
    ProtMatch[i]      = 0.0;
    ProtMatchLevel[i] = 0.0;
    ProtMatchPhase[i] = 0;
  }

  fprintf(stderr,"Reading Blastx data, level... ");
  fflush(stderr);
  
  for( k=0; k<(int)strlen(PAR.getC("blastArg")); k++ ) {
    strcpy(tempname,PAR.getC("fstname"));
    strcat(tempname,".blast");
    i = strlen(tempname);
    tempname[i]   = PAR.getC("blastArg")[k] - 1;
    tempname[i+1] = 0;

    fblast = FileOpen(NULL,tempname, "r");
    if (fblast == NULL) continue;
    
    level      = (PAR.getC("blastArg")[k] - '0') - 1;
    A[0]     = B[0]= 0;
    ProtId   = A;
    PProtId  = B;
    PProtFin = -10000;
    
    while (fscanf(fblast,"%d %d %f %*s %d %s %d %d\n", 
		  &deb, &fin, &score, &phase, ProtId,&ProtDeb,&ProtFin) != EOF) {
      if (abs(fin-deb) > MaxHitLen) {
	fprintf(stderr,"Similarity of extreme length rejected. Check %s\n",ProtId);
	continue;
      }
      
      if (phase < 0) {
	j   = deb;
	deb = fin;
	fin = j;
	j   = ProtDeb;
	ProtDeb = ProtFin;
	ProtFin = j;
      }
      GlobalScore=((REAL)score)/((REAL)abs(fin-deb));
      
      overlap=0;
      // Reconstruction GAPS -> INTRONS
      if ((strcmp(ProtId,PProtId) == 0) && (abs(PProtFin-ProtDeb)<=(MaxOverlap))) {
	// Detection d'un INTRON
	overlap= (PProtFin+1-ProtDeb)*3; // *3 car coord.nucleique
	// overlap >0 : hits chevauchants, overlap <0 : hits espaces
	if ((deb-Pfin+overlap) >= minL8){
	  // Le tableau des introns est rempli prudemment: uniquement les bordures,
	  // et sans serrer pres de l'exon. 
	  for(i = Pfin - (overlap<0) * overlap ;
	      i < Pfin + minL8 - abs(overlap) ; i++) {
	    // debut de l'intron seulement... (dont le score est fonction de
	    // l'exon precedent)
	    if (keyBXLevel[level] >= ProtMatchLevel[i]) {
	      if (keyBXLevel[level] > ProtMatchLevel[i]) {
		ProtMatchLevel[i]= keyBXLevel[level];
		ProtMatch[i]= -PGlobalScore;
		ProtMatchPhase[i]=0;
	      }
	      else {
		if (PGlobalScore >= fabs(ProtMatch[i])){
		  ProtMatch[i]= -PGlobalScore;
		  ProtMatchPhase[i]= 0;
		}
	      }
	    }
	    //		     j=((phase < 0)?-4:4);
	    //		     PlotBarI(i,j,0.6+(level/8.0),1,LevelColor[level]);
	  }
	  for (i = deb - minL8 + abs(overlap) ;
	       i < deb + (overlap<0) * overlap ; i++) {
	    // ...et fin de l'intron (score est fonction de l'exon actuel)
	    if (keyBXLevel[level] >= ProtMatchLevel[i]) {
	      if (keyBXLevel[level] > fabs(ProtMatchLevel[i])) {
		ProtMatchLevel[i]= keyBXLevel[level];
		ProtMatch[i]= -GlobalScore;
		ProtMatchPhase[i]=0;
	      }
	      else
		if (GlobalScore >= fabs(ProtMatch[i])) {
		  ProtMatch[i]= -GlobalScore;
		  ProtMatchPhase[i]= 0;
		}
	    }
	    //		     j=((phase < 0)?-4:4);
	    //		     PlotBarI(i,j,0.6+(level/8.0),1,LevelColor[level]);
	  }
	}
	if (PAR.getI("graph") && level<3) {
	  PlotLine(Pfin,deb,PPhase,phase,0.6+(level/8.0),0.6+(level/8.0),
		   LevelColor[level]);
	  // for(i= Pfin-(overlap<0)*overlap; i < deb+(overlap<0)*overlap ; i++){
	  //   j=((phase < 0)?-4:4);
	  //   PlotBarI(i,j,0.6+(level/8.0),1,LevelColor[level]);
	  // }
	}
      }
      
      PGlobalScore = GlobalScore;
      Pfin = fin;
      tmp  = PProtId;
      PProtId  = ProtId;
      ProtId   = tmp;
      PProtFin = ProtFin;
      PPhase   = phase;
      
      phase = ph06(phase);
      
      // HITS -> CODANT
      if (PAR.getI("graph") && level<3) {
	for (i = deb-1; i < fin; i++) {
	  PlotBarI(i,PhaseAdapt(phase),0.6+(level/8.0),1,LevelColor[level]);
	}
      }
      
      for (i = deb-1; i < fin; i++)
	if (keyBXLevel[level] >= ProtMatchLevel[i])
	  if (keyBXLevel[level] > ProtMatchLevel[i]) {
	    ProtMatchLevel[i] = keyBXLevel[level];
	    ProtMatch[i]= GlobalScore;
	    ProtMatchPhase[i]= PhaseAdapt(phase);
	  }
	  else
	    if (GlobalScore >= fabs(ProtMatch[i])){ 
	      ProtMatch[i] = GlobalScore;
	      ProtMatchPhase[i]= PhaseAdapt(phase);
	    }
    }
    fprintf(stderr,"%d ",level);
    fflush(stderr);
    
    fclose(fblast);       
  }
  fprintf(stderr," done\n");
}

// --------------------------
//  GiveInfo Content BlastX.
// --------------------------
void SensorBlastX :: GiveInfo(DNASeq *X, int pos, DATA *d)
{
  int i;
  for(i=0; i<6; i++)       //exons
    if(ProtMatch[pos] < 0 || 
       ((ProtMatch[pos] > 0) && (ProtMatchPhase[pos] != PhaseAdapt(i))))
      d->ContentScore[i] += -fabs(ProtMatch[pos])*ProtMatchLevel[pos];

  for(i=8; i<13; i++)      //inter & UTRs
    if(ProtMatch[pos] != 0)
      d->ContentScore[i] += -fabs(ProtMatch[pos])*ProtMatchLevel[pos];

  for(i=6; i<8; i++)       //introns
    if(ProtMatch[pos] > 0)
      d->ContentScore[i] += -ProtMatch[pos]*ProtMatchLevel[pos];
}

// -----------------------------------------------
//  Convertit les phases 0-6 en 1 2 3 -1 -2 -3 0.
// -----------------------------------------------
int SensorBlastX :: PhaseAdapt(char p)
{
  if (p >= 12) return 0;
  else if (p < 3) return (1+p);
  else if (p < 6) return (2-p);
  else if (p < 9) return (p-2);
  else return (5-p);
}

// -----------------------------------------------
//  Convertit les phases 1 2 3 -1 -2 -3 0 en 0-6.
// -----------------------------------------------
char SensorBlastX :: ph06(char p)
{
  if (p == 0) return 6;
  else if (p > 0) return (p-1);
  else return 2-p;   
}
