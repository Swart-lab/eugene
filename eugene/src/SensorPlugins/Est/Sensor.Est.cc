#include "Sensor.Est.h"
#include "../MSensor.h"

/*************************************************************
 **                        SensorEst                        **
 *************************************************************/

#define Inconsistent(x) (((x) & Hit) && ((x) & Gap))

extern MasterSensor MS;
extern Parameters   PAR;

int HitsCompare(const void *A, const void *B)
{
  Hits **UA,**UB;
  
  UA = (Hits **) A;
  UB = (Hits **) B;

  if ((*UA)->NGaps > (*UB)->NGaps)   return -1;
  if ((*UA)->NGaps < (*UB)->NGaps)   return  1;
  if ((*UA)->Length > (*UB)->Length) return -1;
  if ((*UA)->Length < (*UB)->Length) return  1;
  return 0;
}

// ----------------------
//  Default constructor.
// ----------------------
SensorEst :: SensorEst (int n) : Sensor(n)
{
  ESTMatch = NULL;
  HitTable = NULL;

  estP = PAR.getD("Est.estP");
  estM = PAR.getD("Est.estM");
}

// ----------------------
//  Default destructor.
// ----------------------
SensorEst :: ~SensorEst ()
{
  if(ESTMatch != NULL)
    delete [] ESTMatch;
  
  if(HitTable != NULL) {
    delete HitTable[NumEST];
    delete HitTable;
  }
  HitTable = NULL;
}

// --------------------------------
//  Init est.
//  Exploiting spliced alignements
//  against EST and complete cDNA.
// --------------------------------
void SensorEst :: Init (DNASeq *X)
{
  FILE *fEST;
  char tempname[FILENAME_MAX+1];
  int  i;
 
  type = Type_Content;
  
  if(ESTMatch != NULL)
    delete [] ESTMatch;
  
  if(HitTable != NULL) {
    delete HitTable[NumEST];
    delete HitTable;
  }
  HitTable = NULL;
 
  ESTMatch = new unsigned char[X->SeqLen+1];
  for (i = 0; i <= X->SeqLen; i++)
    ESTMatch[i] = 0;
 
  strcpy(tempname, PAR.getC("fstname"));
  strcat(tempname, ".est");
  fEST = FileOpen(NULL, tempname, "r");
  NumEST = 0;
  HitTable = ESTAnalyzer(fEST, ESTMatch, estM, &NumEST, X);
  fclose(fEST);
}

// -----------------------
//  GiveInfo signal est.
// -----------------------
void SensorEst :: GiveInfo (DNASeq *X, int pos, DATA *d)
{
  for(int i=0; i<3; i++)              // Exon F
    // Si on a un Gap EST ou si l'on connait le sens du match EST
    if ((ESTMatch[pos] & Gap) || 
	((ESTMatch[pos] & Hit) && !(ESTMatch[pos] & HitForward)))
      d->ContentScore[i] += estP;

  for(int i=3; i<6; i++)              // Exon R
    // Si on a un Gap EST ou si l'on connait le sens du match EST
    if ((ESTMatch[pos] & Gap) ||
	((ESTMatch[pos] & Hit) && !(ESTMatch[pos] & HitReverse)))
      d->ContentScore[i] += estP;

  // Si on a un Hit EST ou si l'on connait le sens du match EST
  if((ESTMatch[pos] & Hit) ||
     ((ESTMatch[pos] & Gap) && !(ESTMatch[pos] & GapForward)))
    d->ContentScore[6] += estP;       // IntronF

  // Si on a un Hit EST ou si l'on connait le sens du match EST
  if((ESTMatch[pos] & Hit) ||
     ((ESTMatch[pos] & Gap) && !(ESTMatch[pos] & GapReverse)))
    d->ContentScore[7] += estP;       // IntronR

  d->ContentScore[8] += ((ESTMatch[pos] & (Gap|Hit)) != 0)*estP;  //InterG

  d->ESTMATCH_TMP = ESTMatch[pos];  // WARNING : EST -> on est dans intron
}

// -----------------------
//  ESTAnalyzer.
// -----------------------
Hits** SensorEst :: ESTAnalyzer(FILE *ESTFile, unsigned char *ESTMatch,
				int EstM, int *NumEST, DNASeq *X)
{
  int i,j,k;
  int Rejected = 0;
  Hits *ThisEST = NULL, *AllEST = NULL;
  Block *ThisBlock = NULL;
  
  fprintf(stderr,"Reading cDNA hits...");
  fflush(stderr);

  AllEST = AllEST->ReadFromFile(ESTFile,NumEST);
  
  fprintf(stderr,"%d sequences read.\n",*NumEST);
  fflush(stderr);

  // on trie les hits sur le nombre de gaps et la
  // longueur. L'idee est d'eliminer les epissages partiels et de
  // favoriser la longueur de toute facon.
  
  Hits **HitTable = new Hits *[*NumEST+1];
  
  for (i = 0, ThisEST = AllEST; i < *NumEST; i++, ThisEST = ThisEST->Next)
    HitTable[i] = ThisEST;

  // pour memorise le premier (a liberer)
  HitTable[*NumEST] = AllEST;
  
  qsort( (void *)HitTable, *NumEST, sizeof(void *), HitsCompare );

  for (int index = 0; index < *NumEST; index++) {
    int Inc;
    int ExonInc,TheStrand;
    REAL WorstSpliceF, WorstSpliceR;
    REAL DonF,AccF,DonR,AccR;
    DATA dTmp;
    
    ThisEST = HitTable[index];

    // Le veritable  brin est a priori indetermine
    TheStrand = HitForward | HitReverse;
    Inc = 0;
    ExonInc = 0;
    WorstSpliceF = WorstSpliceR = 1.0;
        
    // Look for each match in the current Hit
    ThisBlock = ThisEST->Match;

    // First Step: tries to determine strand
    while (ThisBlock) {
      // si ona un gap ?
      if ((ThisBlock->Prev != NULL) && 
	  abs(ThisBlock->Prev->LEnd - ThisBlock->LStart) <= 6) {
	DonF = 0.0;
	DonR = 0.0;
	AccF = 0.0;
	AccR = 0.0;
	
	for (j = -EstM; j <= EstM; j++) {
	  k = Min(X->SeqLen,Max(0,ThisBlock->Prev->End+j+1));
	  MS.GetInfoSpAt(Type_Splice, X, k, &dTmp);
	  DonF = Max(DonF, dTmp.Don[0]);
	  AccR = Max(AccR, dTmp.Acc[1]);
	  
	  k = Min(X->SeqLen,Max(0,ThisBlock->Start+j));
	  if(MS.GetInfoSpAt(Type_Splice, X, k, &dTmp)) {
	    DonR = Max(DonR, dTmp.Don[1]);
	    AccF = Max(AccF, dTmp.Acc[0]);
	  }
	  else {
	    printf("   WARNING: cDNA hits ignored."
		   " No splices sites predicted !\n");
	    return NULL;
	  }
	}

	//	printf("Extreme splices: %f %f - %f %f\n",DonF,AccF,DonR,AccR);
	WorstSpliceF = Min(WorstSpliceF,DonF);
	WorstSpliceF = Min(WorstSpliceF,AccF);
	WorstSpliceR = Min(WorstSpliceR,DonR);
	WorstSpliceR = Min(WorstSpliceR,AccR);
      }
      
      ThisBlock = ThisBlock->Next;
      }
    
    // Tous les blocs ont ete traites
    if (WorstSpliceF == 0.0) TheStrand &= (~HitForward);
    if (WorstSpliceR == 0.0) TheStrand &= (~HitReverse);

    //    printf("Strand %d\n",TheStrand);

    // next iteration on the same EST
    ThisBlock = ThisEST->Match;
    while (TheStrand && ThisBlock) {
      // Check for consistency with already read Hits
      // The Inc flag will keep the Inconsistency status
      // 1 - inconsistent with a previous EST
      // 2 - no splice site on the borders of a gap
      // 3 - an exon contains STRONG donor on both strands

      for (i = ThisBlock->Start+EstM; !Inc && i <= ThisBlock->End-EstM; i++) {

	if (((ESTMatch[i] & Hit) && !(ESTMatch[i] & TheStrand)) |
	    (Inconsistent(ESTMatch[i] | TheStrand))) {
	  fprintf(stderr,"   [%s]: inconsistent hit [%d-%d]\n",
		  ThisEST->Name,ThisBlock->Start+1,ThisBlock->End+1);
	  Inc = 1;
	}
      }
      
      // si on a un gap
      if ((ThisBlock->Prev != NULL) &&
	  abs(ThisBlock->Prev->LEnd - ThisBlock->LStart) <= 6) {
        for (i=ThisBlock->Prev->End+1+EstM; !Inc && i<ThisBlock->Start-EstM; i++) 
	  if (((ESTMatch[i] & Gap) && !(ESTMatch[i] & (TheStrand << HitToGap))) |
	    (Inconsistent(ESTMatch[i] | (TheStrand << HitToGap)))) {
            fprintf(stderr,"   [%s]: inconsistent gap [%d-%d]\n",
                    ThisEST->Name,ThisBlock->Prev->End+2,ThisBlock->Start);
            Inc = 1;
          }
      }

      DonF = 0.0;
      DonR = 0.0;
      
      // calcul des sites d'epissage internes a l'exon
      for (i = ThisBlock->Start+EstM+1; !Inc && i <= ThisBlock->End-EstM-1; i++) {
	MS.GetInfoSpAt(Type_Splice, X, i, &dTmp);
	DonF = Max(DonF, dTmp.Don[0]);
	DonR = Max(DonR, dTmp.Don[1]);
	//	Don[0][i] *= 0.99;
	//	Don[1][i] *= 0.99;
      }
      //      printf("Internal splices: %f - %f\n",DonF,DonR);
      if (DonF > 0.95) ExonInc |= 1;
      if (DonR > 0.95) ExonInc |= 2;
      if (ExonInc == 3 && !Inc && !ThisEST->NGaps) {
	fprintf(stderr,"   [%s]: Gapless EST with strong donor [%d-%d]\n",
		    ThisEST->Name,ThisBlock->Start+1,ThisBlock->End+1);
	Inc = 3;
      }
      
      ThisBlock = ThisBlock->Next;
    }
    
    if (!TheStrand) {
      fprintf(stderr, "   [%s]: no matching splice site\n",ThisEST->Name);
      Inc =2;
    }
    
    // Si une incoherence est detectee, on va jeter la sequence
    
    if (Inc) {
      Rejected++;
      ThisEST->Rejected = TRUE;
      
      ThisBlock = ThisEST->Match;
      while (ThisBlock) {
	for (i = ThisBlock->Start; i<= ThisBlock->End; i++) {
	  if (TheStrand & HitForward) PlotBarI(i, 4,0.9,1,7);
	  if (TheStrand & HitReverse) PlotBarI(i,-4,0.9,1,7);
	}
	
	if (PAR.getI("Output.graph") && (ThisBlock->Prev != NULL) && 
	    abs(ThisBlock->Prev->LEnd-ThisBlock->LStart) <= 6) {
	  if (TheStrand & HitForward) 
	    PlotLine(ThisBlock->Prev->End,ThisBlock->Start,4,4,0.9,0.9,7);
	  if (TheStrand & HitReverse) 
	    PlotLine(ThisBlock->Prev->End,ThisBlock->Start,-4,-4,0.9,0.9,7);
	}
	ThisBlock = ThisBlock->Next;
      }
    }
    // sinon on l'exploite
    else {
      int LBoundary;
      int RBoundary;
      
      ThisBlock = ThisEST->Match;
      
      while (ThisBlock) {
	
	// Aligners tend to extend beyond the true hit on
	// extremities: we remove EstM on frontiers
	
	LBoundary = ((ThisBlock->Prev == NULL) ? 
		     Min(X->SeqLen, ThisBlock->Start+EstM) :
		     ThisBlock->Start);
	
	RBoundary = ((ThisBlock->Next == NULL) ? 
		     Max(0,ThisBlock->End-EstM) :
		     ThisBlock->End);
	
	for (i = LBoundary; i <= RBoundary; i++) {
	  if (ESTMatch[i] & Hit)
	    ESTMatch[i] |= (ESTMatch[i] & TheStrand);
	  else
	    ESTMatch[i] |= TheStrand;
	}

	if (ThisBlock->Prev == NULL)
	  for (i = ThisBlock->Start; i<LBoundary; i++)
	    ESTMatch[i] |= TheStrand << HitToMLeft;

	if (ThisBlock->Next == NULL)
	  for (i = RBoundary+1; i <= ThisBlock->End; i++)
	    ESTMatch[i] |= TheStrand << HitToMRight;
	
	if ((ThisBlock->Prev != NULL) && 
	    abs(ThisBlock->Prev->LEnd-ThisBlock->LStart) <= 6) {
	  
	  for (i = Max(0,ThisBlock->Prev->End+1-EstM); 
	       i < Min(X->SeqLen, ThisBlock->Prev->End+1+EstM); i++) 
	    if (ESTMatch[i] & MLeft)
	      ESTMatch[i] |= (ESTMatch[i] & (TheStrand << HitToMLeft));
	    else
	      ESTMatch[i] |= TheStrand << HitToMLeft;
	  
	  for (i = ThisBlock->Prev->End+1; i < ThisBlock->Start-1; i++)
	    if (ESTMatch[i] & Gap)
	      ESTMatch[i] |= (ESTMatch[i] & (TheStrand << HitToGap));
	    else
	    ESTMatch[i] |= TheStrand << HitToGap;
	  
	  for (i =  Max(0,ThisBlock->Start-1-EstM);
	       i <  Min(X->SeqLen, ThisBlock->Start-1+EstM); i++)
	    if (ESTMatch[i] & MRight)
	      ESTMatch[i] |= (ESTMatch[i] & (TheStrand << HitToMRight));
	    else
	      ESTMatch[i] |= TheStrand << HitToMRight;

	  
	  if (PAR.getI("Output.graph")) {
	    if (TheStrand & HitForward) 
	      PlotLine(ThisBlock->Prev->End,ThisBlock->Start,4,4,0.6,0.6,2);
	    if (TheStrand & HitReverse) 
	      PlotLine(ThisBlock->Prev->End,ThisBlock->Start,-4,-4,0.6,0.6,2);
	  }
	  
	}
	ThisBlock = ThisBlock->Next;
      }
    }
  }
  
  if (PAR.getI("Output.graph"))
    for (i = 0; i < X->SeqLen; i++) {
      if (ESTMatch[i] & HitForward) 
	PlotBarI(i, 4,0.6,1,2);
      if (ESTMatch[i] & HitReverse) 
	PlotBarI(i, -4,0.6,1,2);
    }

  if (Rejected)
    fprintf(stderr,"A total of %d/%d sequences rejected\n",Rejected,*NumEST);
  return HitTable;
}
