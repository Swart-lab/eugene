#include "Sensor.BlastX.h"

/******************************************************************
 **                          SensorBlastX                        **
 ******************************************************************
 * Blastx against protein databases, we simply *
 * enhance coding probability according to phase Dangerous (NR    *
 * contains translation of frameshifted sequences) another        *
 * possibility would be too forbid introns/intergenic states      *
 * 10 levels of confidence may be used.                            *
 * TODO : distinction between intron forward and reverse
 ******************************************************************/

extern Parameters PAR;

// ----------------------
// Default constructor.
// ----------------------
SensorBlastX :: SensorBlastX (int n) : Sensor(n)
{
  keyBXLevel[0] = PAR.getD("BlastX.level0",n);
  keyBXLevel[1] = PAR.getD("BlastX.level1",n);
  keyBXLevel[2] = PAR.getD("BlastX.level2",n);
  keyBXLevel[3] = PAR.getD("BlastX.level3",n);
  keyBXLevel[4] = PAR.getD("BlastX.level4",n);
  keyBXLevel[5] = PAR.getD("BlastX.level5",n);
  keyBXLevel[6] = PAR.getD("BlastX.level6",n);
  keyBXLevel[7] = PAR.getD("BlastX.level7",n);
  keyBXLevel[8] = PAR.getD("BlastX.level8",n);
  keyBXLevel[9] = PAR.getD("BlastX.level9",n);
  blastxM = PAR.getI("BlastX.blastxM",n);
  minIn = PAR.getI("EuGene.minIn");
}

// ----------------------
//  Default destructor.
// ----------------------
SensorBlastX :: ~SensorBlastX ()
{
  vPos.clear();
  vPMatch.clear();
  vPMLevel.clear();
  vPMPhase.clear();
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
  int   deb, fin, phase, ProtDeb, ProtFin, PProtFin, level;
  float score;
  int   Pfin = 0, Pphase = 0;
  char  A[128], B[128];
  char *ProtId, *PProtId, *tmp;
  REAL GlobalScore;
  REAL PGlobalScore = 0;
  const int MaxHitLen  = 15000;
  char tempname[FILENAME_MAX+1];

  type = Type_Content;

  index = 0;

  vPos.clear();
  vPMatch.clear();
  vPMLevel.clear();
  vPMPhase.clear();

  ProtMatch      = new REAL[Len+1];
  ProtMatchLevel = new REAL[Len+1];
  ProtMatchPhase = new int[Len+1];
  for (i = 0; i<= Len; i++) {
    ProtMatch[i]      = 0.0;
    ProtMatchLevel[i] = 0.0;
    ProtMatchPhase[i] = 0;
  }

  fprintf(stderr,"Reading Blastx data, level...");
  fflush(stderr);

  for( k=0; k<(int)strlen(PAR.getC("BlastX.levels",GetNumber())); k++ ) {
  // for each level given in arg: 
  // e.g. with -b09 : 1st: level 0 (k=0), 2nd: level 9 (k=1)
    strcpy(tempname,PAR.getC("fstname"));
    strcat(tempname,".blast");
    i = strlen(tempname);
    tempname[i]   = PAR.getC("BlastX.levels",GetNumber())[k];
    tempname[i+1] = 0;
    // check the corresponding .blastN file (N=level given in arg)
    fblast = FileOpen(NULL,tempname, "r");
    if (fblast == NULL) {
      fprintf(stderr, "BlastX error : unable to open file %s\n",tempname);
      exit(1);
    }

    level      = (PAR.getC("BlastX.levels",GetNumber())[k] - '0');
    A[0]     = B[0]= 0;
    ProtId   = A;
    PProtId  = B;
    PProtFin = -10000;

    // Read the blast file
    while (fscanf(fblast,"%d %d %f %*s %d %s %d %d\n", 
		  &deb, &fin, &score, &phase, ProtId,&ProtDeb,&ProtFin) != EOF) {
      // for each HSP
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

      // GAPS -> INTRONS ; the "intron" consideration between 2 HSP requires:
      // 1: same prot ID, 2: close prot boundaries (blastxM param), 3: same nucleotidic strand...
      if ( (strcmp(ProtId,PProtId) == 0) && (abs(PProtFin-ProtDeb)<=(blastxM)) && (phase*Pphase>0) ) {
	// INTRON
	if ((deb-Pfin) >= minIn) {
	  for(i = Pfin;  i < Pfin+(minIn)/2;i++) {
	    // begining of the intron only, during a short region
	    // the score depends on the adjacent exon
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
	  }
	  for (i = deb - minIn/2 ; i< deb ; i++) {
	    // ...and the end of the intron
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
	if (PAR.getI("Output.graph") && k <3) {
	  PlotLine(Pfin,deb,Pphase,phase,0.6+(k/8.0),0.6+(k/8.0),
		   LevelColor[k]);
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
      Pphase   = phase;

      phase = ph06(phase);

      // HITS -> CODING
      if (PAR.getI("Output.graph") && k<3) {
	for (i = deb-1; i < fin; i++) {
	  PlotBarI(i,PhaseAdapt(phase),0.6+(k/8.0),1,LevelColor[k]);
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

  for (i = 0; i<= Len; i++)
    if(ProtMatch[i] != 0.0) {
      vPos.push_back    ( i );
      vPMatch.push_back ( ProtMatch[i] );
      vPMLevel.push_back( ProtMatchLevel[i] );
      vPMPhase.push_back( ProtMatchPhase[i] );
    }

  delete [] ProtMatch;
  delete [] ProtMatchLevel;
  delete [] ProtMatchPhase;

  fprintf(stderr,"done\n");
}

// --------------------------
//  GiveInfo Content BlastX.
// --------------------------
void SensorBlastX :: GiveInfo(DNASeq *X, int pos, DATA *d)
{
  int i;

  if((index != 0                &&  vPos[index-1] >= pos) ||
     (index < (int)vPos.size()  &&  vPos[index]   <  pos))
    {
      iter = lower_bound(vPos.begin(), vPos.end(), pos);
      if(*iter == pos) {
	for(i=0; i<6; i++)       //exons
	  if(vPMatch[iter-vPos.begin()] < 0 ||
	     ((vPMatch[iter-vPos.begin()] > 0) && (vPMPhase[iter-vPos.begin()] != PhaseAdapt(i))))
	    d->contents[i] +=  -fabs(vPMatch[iter-vPos.begin()])*vPMLevel[iter-vPos.begin()];
	
	for(i=8; i<13; i++)      //inter & UTRs
	  if(vPMatch[iter-vPos.begin()] != 0)
	    d->contents[i] +=  -fabs(vPMatch[iter-vPos.begin()])*vPMLevel[iter-vPos.begin()];
	
	for(i=6; i<8; i++)       //introns
	  if(vPMatch[iter-vPos.begin()] > 0)
	    d->contents[i] +=  -vPMatch[iter-vPos.begin()]*vPMLevel[iter-vPos.begin()];
	index = iter-vPos.begin() + 1;
      }
      else index = iter-vPos.begin();
    }
  else if( index < (int)vPos.size()  &&  vPos[index] == pos )
    {
      for(i=0; i<6; i++)       //exons
	if(vPMatch[index] < 0 || 
	   ((vPMatch[index] > 0) && (vPMPhase[index] != PhaseAdapt(i))))
	  d->contents[i] +=  -fabs(vPMatch[index])*vPMLevel[index];
      
      for(i=8; i<13; i++)      //inter & UTRs
	if(vPMatch[index] != 0)
	  d->contents[i] +=  -fabs(vPMatch[index])*vPMLevel[index];
      
      for(i=6; i<8; i++)       //introns
	if(vPMatch[index] > 0)
	  d->contents[i] += -vPMatch[index]*vPMLevel[index];
      index++;
    }
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

// ----------------------------
//  Plot Sensor information
// ----------------------------
void SensorBlastX :: Plot(DNASeq *X)
{
}

// ------------------
//  Post analyse
// ------------------
void SensorBlastX :: PostAnalyse(Prediction *pred)
{
  int state= 0, stateNext = 0, stateBack = 0;
  int pos= 0,   posNext   = 0, posBack   = 0;
  int begCDS= 0, endCDS =0;
  int CodingNuc = 0;
  int SupportedNuc = 0;

  if (!PAR.getI("BlastX.PostProcess",GetNumber())) return;

  index=0;
  for(int i=pred->size()-1; i!=-1; i--) {
    if(i != pred->size()-1) {
      stateBack = pred->getState(i+1);
      posBack   = pred->getPos(i+1);
    }
    state = pred->getState(i);
    pos   = pred->getPos(i);
    if(i != 0) {
      stateNext = pred->getState(i-1);
      posNext   = pred->getPos(i-1);
    }

    if (state <= ExonR3) { // EXON

      if (stateBack >= UTR5F) { // c'est le premier (ou dernier) d'une CDS
	CodingNuc = 0;
	SupportedNuc = 0;
	begCDS=posBack+1;
      }

      CodingNuc += pos-posBack;

      while ( (index < (int)vPos.size()) && (vPos[index] < posBack)) index++;
      // index=indice du premier match etant > ou = au debut de l'exon
      // fprintf(stderr,"premier match > ou = au debut de l'exon: pos %d, index=%d\n",vPos[index],index);

      if (index < (int)vPos.size()) { // il reste des hits
	// pour chaque nuc de l'exon supporte par un hit
	while (vPos[index]<=pos  &&  index < (int)vPos.size()) {
	  if (State2Phase[state] ==  vPMPhase[index]) SupportedNuc++;
	  index++;
	}
      }
      
      if ( ((state <= ExonR3) && (stateNext >= InterGen5)) || // fin de gene
	   (i==0 && state <= IntronR3)) {                 // ou fin seq(gene en cours)
	endCDS = pos;
	printf("      CDS (prot.)  %7d %7d    %5d     supported on %d/%d coding.\n",begCDS,endCDS,endCDS-begCDS+1,SupportedNuc,CodingNuc);
      }
    }
  }
}
