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
  int i,k;
  char   tempname[FILENAME_MAX+1];
  FILE * fblast;

  N        = n;
  type     = Type_Content;
  HitTable = NULL;
  ppNumber = PAR.getI("BlastX.PPNumber",N);
  stepid   = PAR.getI("Output.stepid");
  minIn    = PAR.getI("BlastX.minIn");
  levels   = PAR.getC("BlastX.levels",  N);  
 
  Hits   *AllProt = NULL;
  NumProt = 0;
  
  fprintf(stderr,"Reading BlastX data, level...");
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
    
    AllProt = AllProt->ReadFromFile(fblast, &NumProt, (levels[k] - '0'), 20);
    fprintf(stderr,"%c ",levels[k]);
    fflush(stderr);
    
    fclose(fblast);       
  }
  fprintf(stderr,"done\n");
  
  
  Hits *ThisProt = AllProt;
  HitTable = new Hits *[NumProt+1];
  for (i=0;  i<NumProt;  i++, ThisProt=ThisProt->Next)
    HitTable[i] = ThisProt;
  
  //for (i=0;  i<NumProt;  i++)
  //printf("Name:%s\tLevel:%d\n", HitTable[i]->Name,HitTable[i]->Level);
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
  int i, j;
  int Len = X->SeqLen;
  const int LevelColor[3] = {6,7,8}; 
  int   level, levelidx, Pphase = 0;
  double GlobalScore,    PGlobalScore = 0;
 
  vPos.clear();
  vPMatch.clear();
  vPMLevel.clear();
  vPMPhase.clear();

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
  ProtMatch      = new double[Len+1];
  ProtMatchLevel = new double[Len+1];
  ProtMatchPhase = new int[Len+1];

  for (i=0; i<= Len; i++) {
    ProtMatch[i]      = 0.0;
    ProtMatchLevel[i] = 0.0;
    ProtMatchPhase[i] = 0;
  }

  for (i=0;  i<NumProt;  i++) {

    Block *MyHSP = HitTable[i]->Match;
    level = HitTable[i]->Level;
    levelidx = rindex(levels,'0'+level)-levels;

    while (MyHSP != NULL) {
      
      GlobalScore = (double)(MyHSP->Score) / (double)abs(MyHSP->End-MyHSP->Start);
      
      // GAPS -> INTRONS ; the "intron" consideration between 2 HSP
      // requires close prot boundaries (blastxM param)
      if ( MyHSP->Prev &&
	   (abs(MyHSP->Prev->LEnd-MyHSP->LStart) <= blastxM )) {
	
	// INTRON
	if ((MyHSP->Start-MyHSP->Prev->End) >= minIn) {

	  for(j = MyHSP->Prev->End+1; j < MyHSP->Prev->End+1+(minIn)/2; j++) {
	    // begining of the intron only, during a short region
	    // the score depends on the adjacent exon
	    if (keyBXLevel[level] >= ProtMatchLevel[j]) {
	      if (keyBXLevel[level] > ProtMatchLevel[j]) {
		ProtMatchLevel[j]= keyBXLevel[level];
		ProtMatch[j]= -PGlobalScore;
		ProtMatchPhase[j]=0;
	      }
	      else {
		if (PGlobalScore >= fabs(ProtMatch[j])){
		  ProtMatch[j]= -PGlobalScore;
		  ProtMatchPhase[j]= 0;
		}
	      }
	    }
	  }

	  for (j = MyHSP->Start+1 - minIn/2; j < MyHSP->Start+1; j++) {
	    // ...and the end of the intron
	    if (keyBXLevel[level] >= ProtMatchLevel[j]) {
	      if (keyBXLevel[level] > fabs(ProtMatchLevel[j])) {
		ProtMatchLevel[j]= keyBXLevel[level];
		ProtMatch[j]= -GlobalScore;
		ProtMatchPhase[j]=0;
	      }
	      else
		if (GlobalScore >= fabs(ProtMatch[j])) {
		  ProtMatch[j]= -GlobalScore;
		  ProtMatchPhase[j]= 0;
		}
	    }
	    //	;
	    //	PlotBarI(i,((phase < 0)?-4:4),0.6+(level/8.0),1,LevelColor[level]);
	  }
	}

	if (PAR.getI("Output.graph") && levelidx <3) {
	  PlotLine(MyHSP->Prev->End,MyHSP->Start,Pphase,MyHSP->Phase,
		   0.6+(levelidx/8.0),0.6+(levelidx/8.0),LevelColor[levelidx]);
	  //for(i=Pfin-(overlap<0)*overlap; i < deb+(overlap<0)*overlap ; i++){
	  //   j=((phase < 0)?-4:4);
	  //   PlotBarI(i,j,0.6+(level/8.0),1,LevelColor[level]);
	  //}
	}
      }

      // HITS -> CODING
      if (PAR.getI("Output.graph") && levelidx<3) {
	for (j = MyHSP->Start-1; j < MyHSP->End; j++) {
	  PlotBarI(j,MyHSP->Phase,
		   0.6+(levelidx/8.0),1,LevelColor[levelidx]);
	}
      }
      
      for (j = MyHSP->Start; j < MyHSP->End+1; j++)
	if (keyBXLevel[level] >= ProtMatchLevel[j])
	  if (keyBXLevel[level] > ProtMatchLevel[j]) {
	    ProtMatchLevel[j] = keyBXLevel[level];
	    ProtMatch[j]= GlobalScore;
	    ProtMatchPhase[j]= MyHSP->Phase;
	  }
	  else
	    if (GlobalScore >= fabs(ProtMatch[j])){ 
	      ProtMatch[j] = GlobalScore;
	      ProtMatchPhase[j]= MyHSP->Phase;
	    }


      PGlobalScore = GlobalScore;
      Pphase = MyHSP->Phase;
      MyHSP = MyHSP->Next;
    }
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

  index = 0;
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
      
      if (pprocess == 1) {
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
      }

      // Si fin de gene ou fin seq gene en cours
      if ((state <= TermR3 && stateNext >= InterGen) ||
	  (i==0 && (state <= IntronR3))) {
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
    end   = pred->getPos(i);
    if(i != pred->size()-1) start = pred->getPos(i+1) + 1;
    else                    start = 1;
     
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
		   (((NumGene-1)*stepid)+1),numF,fea,start,end,len,strand);
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
	    for(j=0; j<NumProt && j<ppNumber && TMPHitTable[j]->Support!=0;j++)
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
	     (((NumGene-1)*stepid)+1), fea, start, end, len, strand);
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
      for(j=0; j<NumProt && j<ppNumber && TMPHitTable[j]->Support!=0; j++)   
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
      from = Max(beg, ThisBlock->Start+1);
      to   = Min(end, ThisBlock->End+1);
       
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
  delete [] Sup;
  return supported;
}
