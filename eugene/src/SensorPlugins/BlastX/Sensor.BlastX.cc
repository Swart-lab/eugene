/*****************************************************************************/
/*             Copyright (c) 2004 by INRA. All rights reserved.              */
/*                 Redistribution is not permitted  without                  */
/*                 the express written permission of  INRA.                  */
/*                   Mail : eugene@ossau.toulouse.inra.fr                    */
/*---------------------------------------------------------------------------*/
/* File         : EuGeneTk/SensorPlugins/BlastX/Sensor.BlastX.cc             */
/* Description  : Sensor BlastX                                              */
/* Authors      : P.Bardou, S.Foissac, M.J.Cros, A.Moisan, T.Schiex          */
/* History      : July 2004                                                  */
/*****************************************************************************/

#include "Sensor.BlastX.h"

/******************************************************************
 **                          SensorBlastX                        **
 ******************************************************************
 * Blastx against protein databases, we simply                    *
 * enhance coding probability according to phase Dangerous (NR    *
 * contains translation of frameshifted sequences) another        *
 * possibility would be too forbid introns/intergenic states      *
 * 10 levels of confidence may be used.                           *
 *****************************************************************/

extern Parameters PAR;

int HitsCompareLex(const void *A, const void *B)
{
  Hits **UA,**UB;
  
  UA = (Hits **) A;
  UB = (Hits **) B;
  
  if ((*UA)->Start < (*UB)->Start) return -1;
  if ((*UA)->Start > (*UB)->Start) return 1;
  if ((*UA)->NGaps > (*UB)->NGaps) return -1;
  if ((*UA)->NGaps < (*UB)->NGaps) return 1;
  if ((*UA)->Length > (*UB)->Length) return -1;
  if ((*UA)->Length < (*UB)->Length) return 1;
  return 0;
}

int HitsCompareSup(const void *A, const void *B)
{
  Hits **UA,**UB;
  
  UA = (Hits **) A;
  UB = (Hits **) B;
  
  if ((*UA)->Support > (*UB)->Support) return -1;
  if ((*UA)->Support < (*UB)->Support) return 1;
  return 0;
}

// ----------------------
// Default constructor.
// ----------------------
SensorBlastX :: SensorBlastX (int n, DNASeq *X) : Sensor(n)
{
  type     = Type_Content;
  HitTable = NULL;
  N        = n;
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
  if(HitTable != NULL) delete [] HitTable;
  HitTable = NULL;
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
  char  A[128],   B[128];
  char  *ProtId,  *PProtId, *tmp, *levels;
  double GlobalScore;
  double PGlobalScore = 0;
  const int MaxHitLen = 15000;
  char   tempname[FILENAME_MAX+1];

  // For postprocess 2
  Hits   *AllProt = NULL;
  NumProt = 0;

  vPos.clear();
  vPMatch.clear();
  vPMLevel.clear();
  vPMPhase.clear();

  if(HitTable != NULL) delete [] HitTable;
  HitTable = NULL;

  keyBXLevel[0] = PAR.getD("BlastX.level0*", N);
  keyBXLevel[1] = PAR.getD("BlastX.level1*", N);
  keyBXLevel[2] = PAR.getD("BlastX.level2*", N);
  keyBXLevel[3] = PAR.getD("BlastX.level3*", N);
  keyBXLevel[4] = PAR.getD("BlastX.level4*", N);
  keyBXLevel[5] = PAR.getD("BlastX.level5*", N);
  keyBXLevel[6] = PAR.getD("BlastX.level6*", N);
  keyBXLevel[7] = PAR.getD("BlastX.level7*", N);
  keyBXLevel[8] = PAR.getD("BlastX.level8*", N);
  keyBXLevel[9] = PAR.getD("BlastX.level9*", N);
  blastxM       = PAR.getI("BlastX.blastxM*",N);
  ppNumber      = PAR.getI("BlastX.PPNumber",N);
  minIn         = PAR.getI("BlastX.minIn");
  levels        = PAR.getC("BlastX.levels",  N);
  ProtMatch      = new double[Len+1];
  ProtMatchLevel = new double[Len+1];
  ProtMatchPhase = new int[Len+1];
  for (i=0; i<= Len; i++) {
    ProtMatch[i]      = 0.0;
    ProtMatchLevel[i] = 0.0;
    ProtMatchPhase[i] = 0;
  }

  fprintf(stderr,"Reading Blastx data, level...");
  fflush(stderr);
  for( k=0; k<(int)strlen(levels); k++ ) {
    // for each level given in arg: 
    // e.g. with -b09 : 1st: level 0 (k=0), 2nd: level 9 (k=1)
    strcpy(tempname,PAR.getC("fstname"));
    strcat(tempname,".blast");
    i = strlen(tempname);
    tempname[i]   = levels[k];
    tempname[i+1] = 0;
    // check the corresponding .blastN file (N=level given in arg)
    fblast = FileOpen(NULL,tempname, "r");
    if (fblast == NULL) {
      fprintf(stderr, "BlastX error : unable to open file %s\n",tempname);
      exit(1);
    }
 
    level    = (levels[k] - '0');
    A[0]     = B[0]= 0;
    ProtId   = A;
    PProtId  = B;
    PProtFin = -10000;

    // For postprocess 2
    AllProt = AllProt->ReadFromFile(fblast, &NumProt, level, 20);
    fseek(fblast, 0, SEEK_SET);

    // Read the blast file
    while (fscanf(fblast,"%d %d %f %*s %d %s %d %d\n", &deb, &fin, &score,
		  &phase, ProtId,&ProtDeb,&ProtFin) != EOF) {
      // for each HSP
      if (abs(fin-deb) > MaxHitLen) {
	fprintf(stderr,"Similarity of extreme length rejected. Check %s\n",
		ProtId);
	continue;
      }

      j       = Min(deb,fin);
      fin     = Max(deb,fin);
      deb     = j;
      j       = Min(ProtDeb,ProtFin);
      ProtFin = Max(ProtDeb,ProtFin);
      ProtDeb = j;

      if (phase < 0) {
	j       = ProtDeb;
	ProtDeb = ProtFin;
	ProtFin = j;
      }
      GlobalScore = ((double)score) / ((double)abs(fin-deb));

      // GAPS -> INTRONS ; the "intron" consideration between 2 HSP requires:
      // 1: same prot ID
      // 2: close prot boundaries (blastxM param)
      // 3: same nucleotidic strand...
      if ( (strcmp(ProtId,PProtId) == 0) && 
	   (abs(PProtFin-ProtDeb) <= (blastxM)) && (phase*Pphase > 0) ) {
	// INTRON
	if ((deb-Pfin) >= minIn) {
	  for(i = Pfin; i < Pfin+(minIn)/2; i++) {
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
	  for (i = deb - minIn/2; i<deb; i++) {
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
	    //	j=((phase < 0)?-4:4);
	    //	PlotBarI(i,j,0.6+(level/8.0),1,LevelColor[level]);
	  }
	}
	if (PAR.getI("Output.graph") && k <3) {
	  PlotLine(Pfin,deb,Pphase,phase,0.6+(k/8.0),0.6+(k/8.0),
		   LevelColor[k]);
	  //for(i=Pfin-(overlap<0)*overlap; i < deb+(overlap<0)*overlap ; i++){
	  //   j=((phase < 0)?-4:4);
	  //   PlotBarI(i,j,0.6+(level/8.0),1,LevelColor[level]);
	  //}
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

  index = 0;

  // For postprocess 2
  Hits *ThisProt = AllProt;
  HitTable = new Hits *[NumProt+1];
  for (i=0;  i<NumProt;  i++, ThisProt=ThisProt->Next)
    HitTable[i] = ThisProt;

  //for (i=0;  i<NumProt;  i++)
  //printf("Name:%s\tLevel:%d\n", HitTable[i]->Name,HitTable[i]->Level);
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
	     ((vPMatch[iter-vPos.begin()] > 0) &&
	      (vPMPhase[iter-vPos.begin()] != PhaseAdapt(i))))
	    d->contents[i] += -fabs(vPMatch[iter-vPos.begin()])*vPMLevel[iter-vPos.begin()];
	
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
  int state  = 0, stateNext = 0, stateBack = 0;
  int pos    = 0, posNext   = 0, posBack   = 0;
  int begCDS = 1, endCDS    = 0;
  int CodingNuc    = 0;
  int SupportedNuc = 0;
  int NumGene      = 1;  

  int pprocess = PAR.getI("BlastX.PostProcess",GetNumber());
  if (!pprocess)    return;
  if (NumProt == 0) return;

  index=0;
  qsort((void *)HitTable, NumProt, sizeof(void *), HitsCompareLex);
  // Reset static in Prot Support or FEA Support
  ProtSupport(NULL, 100, 0, NULL, 0, 0);
  
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
    
    if (state <= TermR3) { // EXON
      if (stateBack >= UTR5F) { // c'est le premier (ou dernier) d'une CDS
	CodingNuc = 0;
	SupportedNuc = 0;
	begCDS=posBack+1;
	NumGene++;
      }
      
      CodingNuc += pos-posBack;
      
      while ( (index < (int)vPos.size()) && (vPos[index] < posBack)) index++;
      // index=indice du premier match etant > ou = au debut de l'exon
      // fprintf(stderr,"premier match > ou = au debut de l'exon: pos %d,
      // index=%d\n",vPos[index],index);
      
      if (index < (int)vPos.size()) { // il reste des hits
	// pour chaque nuc de l'exon supporte par un hit
	while (index < (int)vPos.size() && vPos[index]<pos) {
	  if (State2Phase[state] ==  vPMPhase[index]) SupportedNuc++;
	  index++;
	}
      }

      // Si fin de gene ou fin seq gene en cours
      if ((state <= TermR3 && stateNext >= InterGen) ||
	  (i==0 && (state <= IntronR3 || state == UTR3R || state == UTR5F))) {
	endCDS = pos;
	if (pprocess == 1) {
	  printf("      CDS (prot.)  %7d %7d    %5d     ",
		 begCDS,endCDS,endCDS-begCDS+1);
	    printf("supported on %d/%d coding.\n", SupportedNuc,CodingNuc);
	}
	else {
	  ProtSupport(pred, begCDS, endCDS, HitTable, NumProt, NumGene);
	}
      }
    }
  }
}

// -------------------------------------------------------------------------
//  Post analyse (Analyse de la prédiction (features) par rapport aux Prot)
//    debut/fin  = debut/fin de traduit
// -------------------------------------------------------------------------
void SensorBlastX :: ProtSupport(Prediction *pred, int debut, int fin,
				 Hits **HitTable,  int Size,  int NumGene)
{
  static int ProtIndex;
  Block *ThisBlock;
  int SupportProt, i, j;
  int from  = 0, to = 0;
  int state = 0;
  std::vector <int> vSupProtI;                // index prot supportant pred
  Hits **TMPHitTable = new Hits *[NumProt+1]; // need for sort by % support

  if (pred == NULL) {
    ProtIndex = 0;
    return;
  }
  
  /***********************************************************************/
  /** Objectif : obtenir un vecteur contenant les index des prots       **/
  /**            qui supportent la prédiction                           **/
  /***********************************************************************/
  // si l'iteration precedente a atteint l'extremite 
  if (ProtIndex >= Size)  ProtIndex = Max(0, Size-1);

  // on rembobine....
  while ((ProtIndex > 0) && (HitTable[ProtIndex]->End > debut))
    ProtIndex--;

  if (ProtIndex >= 0  &&  HitTable[ProtIndex]->End < debut) ProtIndex++;
  
  while (ProtIndex >=0  &&  ProtIndex < Size) {
    // la derniere prot exploitable est passee
    if (HitTable[ProtIndex]->Start > fin) break;

    SupportProt = 0;
    ThisBlock = HitTable[ProtIndex]->Match;
    
    while (ThisBlock) {
      from = Max(debut, ThisBlock->Start);
      to   = Min(fin,   ThisBlock->End);
      if (from < to)
	SupportProt = 1;
      ThisBlock = ThisBlock->Next;
    }
    
    // Si Prot supporte
    if (SupportProt)
      vSupProtI.push_back( ProtIndex );
    ProtIndex++;
  }
  
  
  /***********************************************************************/
  /* Objectif : Analyser chaque feature predite -> supportée ?           */
  /***********************************************************************/
  int  start = 0, end = 0, len;
  char fea[5];
  char strand = '+';
  if (PhaseAdapt(pred->getStateForPos(debut+1)) < 0) strand = '-';
  int NumFEA    = 0;
  int codingNuc = 0;
  if(strand == '-') NumFEA = pred->nbExon(NumGene) + 1;
    
  for(i=pred->size()-1; i!=-1; i--) {
    state = pred->getState(i);
    start = pred->getPos(i+1) + 1;
    end   = pred->getPos(i);
    if(end > debut  &&  end <= fin+1)
      {
	len      = 0;
	int numF = -1;

	if(state <= TermR3) {
	  (strand == '-') ? NumFEA-- : NumFEA++;
	  strcpy(fea, "Exon");
	  numF = NumFEA;
	  codingNuc += end - start + 1;
	}
	else if(state == UTR5F  ||  state == UTR5R  ||
		state == UTR3F  ||  state == UTR3R) {
	  strcpy(fea, "UTR");
	  numF      = 0;
	}

      	if(numF != -1) {
	  if((int)vSupProtI.size() > 0)
	    // Longueur totale supportée par les prots
	    len = LenSup(HitTable, pred, vSupProtI, -1, start, end);
	  
	  if(len > 0) {
	    printf("%d.%d\tEuGene_prot\t%s\t%d\t%d\t%d\t%c\t.\t",
		   NumGene, numF, fea, start, end, len, strand);
	    printf("%d\t", (int)((float)len/(end-start+1)*100));
	    
	    for (int k=0;  k<NumProt;  k++)
	      HitTable[k]->Support = 0;

	    for(j=0; j<(int)vSupProtI.size(); j++) {
	      // Longueur supportée par la prot
	      len = LenSup(HitTable, pred, vSupProtI, j, start, end);
	      HitTable[vSupProtI[j]]->Support = (int)((float)len/(end-start+1)*100);
	    }
	    // On copie la hittable pour trier sur le % supporté
	    for (int k=0; k<NumProt; k++)
	      TMPHitTable[k] = HitTable[k];
	    qsort((void*)TMPHitTable, NumProt, sizeof(void*), HitsCompareSup);
	   
	    // On affiche les ppNumber premiers hits supportant
	    for(j=0; j<ppNumber && TMPHitTable[j]->Support!=0; j++)
	      printf("%s(%d,%d) ", TMPHitTable[j]->Name,
		     TMPHitTable[j]->Support, TMPHitTable[j]->Level);
	    printf("\n");
	  }
	}
      }
  }

  /***********************************************************************/
  /* Objectif : Analyser la CDS et le gene predit -> supportée ?         */
  /***********************************************************************/
  start = debut;
  end = fin;
  strcpy(fea, "CDS");

  if (end >= start) {
    len = 0;
    if((int)vSupProtI.size() > 0)
      // Longueur totale supportée par les prots
      len = LenSup(HitTable, pred, vSupProtI, -1, start, end);
    
    if(len > 0) {
      printf("%d\tEuGene_prot\t%s\t%d\t%d\t%d\t%c\t.\t",
	     NumGene, fea, start, end, len, strand);
      printf("%d\t", (int)((float)len/codingNuc*100));
      for (int k=0;  k<NumProt;  k++)
	HitTable[k]->Support = 0;
	
      for(j=0; j<(int)vSupProtI.size(); j++) {
	// Longueur supportée par la prot
	len = LenSup(HitTable, pred, vSupProtI, j, start, end);
	HitTable[vSupProtI[j]]->Support = (int)((float)len/codingNuc*100);
      }
      // On copie la hittable pour trier sur le % supporté
      for (int k=0; k<NumProt; k++)
	TMPHitTable[k] = HitTable[k];
      qsort((void*)TMPHitTable, NumProt, sizeof(void*), HitsCompareSup);
    
      // On affiche les ppNumber premiers hits supportant
      for(j=0; j<ppNumber && TMPHitTable[j]->Support!=0; j++)   
	printf("%s(%d,%d) ", TMPHitTable[j]->Name,
	       TMPHitTable[j]->Support, TMPHitTable[j]->Level);
      printf("\n");
    }
  }
  vSupProtI.clear();
  if(TMPHitTable != NULL) delete [] TMPHitTable;
  TMPHitTable = NULL;
  return;
}

// -------------------------------------------------------------------------
//  Lenght supported by Prot.
//    If index = -1 -> Lenght supported by ALL Prot.
//    Else          -> Lenght supported by ONE Prot.
// -------------------------------------------------------------------------
int SensorBlastX :: LenSup(Hits **HitTable, Prediction *pred,
			   std::vector<int> vSupProtI,
			   int index, int beg, int end)
{
  int supported = 0;
  unsigned char *Sup;
  Block *ThisBlock;
  int from  = 0, to = 0;
  int state = 0;
  int i, j;

  Sup = new unsigned char[end-beg+1]; 
  for (i=0; i<=end-beg; i++)
    Sup[i]=0;
  
  index == -1  ?  i=0 : i=index;
  for(; i<(int)vSupProtI.size() || i==index; i++) {
    ThisBlock = HitTable[vSupProtI[i]]->Match;
    while (ThisBlock) {
      from = Max(beg, ThisBlock->Start);
      to   = Min(end, ThisBlock->End);
       
      for (j=from; j<=to; j++) {
	state = pred->getStateForPos(j);
	if (!Sup[j-beg]  &&  State2Phase[state] == ThisBlock->Phase) {
	  Sup[j-beg] = TRUE;
	  supported++;
	}
      }
      ThisBlock = ThisBlock->Next;
    }
    if(index != -1) break;
  }
  delete Sup;
  return supported;
}
