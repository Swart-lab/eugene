#include "Sensor.Est.h"
#include "../../EuGene/MSensor.h"

/*************************************************************
 **                        SensorEst                        **
 *************************************************************/
#define EXTREMEMARGIN 1
#define Inconsistent(x) (((x) & Hit) && ((x) & Gap))

extern MasterSensor *MS;
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
  return strcmp((*UA)->Name,(*UB)->Name); //pour un tri sans ambiguité
}

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

// ----------------------
//  Default constructor.
// ----------------------
SensorEst :: SensorEst (int n, DNASeq *X) : Sensor(n)
{
  type = Type_Content;
  HitTable = NULL;
  N = n;
}

// ----------------------
//  Default destructor.
// ----------------------
SensorEst :: ~SensorEst ()
{
  vPos.clear();
  vESTMatch.clear();
  if(HitTable != NULL) delete [] HitTable;
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
 
  vPos.clear();
  vESTMatch.clear();
  
  if(HitTable != NULL) delete [] HitTable;
  HitTable = NULL;

  estP = PAR.getD("Est.estP*",N);
  estM = PAR.getI("Est.estM",N);
  utrP = PAR.getD("Est.utrP*",N);
  utrM = PAR.getI("Est.utrM",N);
  spliceBoost = PAR.getD("Est.SpliceBoost*",N);
  DonorThreshold = PAR.getD("Est.StrongDonor*");
  DonorThreshold = log(DonorThreshold/(1-DonorThreshold));

  index = 0;
 
  ESTMatch = new unsigned char[X->SeqLen+1];
  for (i = 0; i <= X->SeqLen; i++)
    ESTMatch[i] = 0;
 
  strcpy(tempname, PAR.getC("fstname"));
  strcat(tempname, ".est");
  fEST = FileOpen(NULL, tempname, "r");
  NumEST = 0;
  HitTable = ESTAnalyzer(fEST, ESTMatch, estM, &NumEST, X);
  fclose(fEST);
  
  for (i = 0; i<= X->SeqLen; i++)
    if(ESTMatch[i] != 0) {
      vPos.push_back      ( i );
      vESTMatch.push_back ( ESTMatch[i] );
    }
  delete [] ESTMatch;

  //for(int jj=0;jj<(int)vPos.size();jj++)
  //printf("vPos[%d]:%d\tvESTM[%d]:%d\n",jj,vPos[jj]+1,jj,vESTMatch[jj]);
}

// -----------------------
//  GiveInfo signal est.
// -----------------------
void SensorEst :: GiveInfo (DNASeq *X, int pos, DATA *d)
{
  unsigned char cESTMatch = 0; // current ESTMatch
  
  // Peut on faire un bete acces sequentiel ?
  if((index != 0                &&  vPos[index-1] >= pos) ||
     (index < (int)vPos.size()  &&  vPos[index]   <  pos))  {
    // Non... on se repositionne en temps logarithmique
    iter = lower_bound(vPos.begin(), vPos.end(), pos);
    index = iter-vPos.begin();
  }
  // On est juste avant ou sur pos
  
  // Si on est dessus
  if (index < (int)vPos.size()  &&  vPos[index] == pos) {
    
    cESTMatch = vESTMatch[index];

    // Favor splice sites in marginal exon-intron regions
    if ((cESTMatch & MLeftForward) && 
	d->sig[DATA::Don].weight[Signal::Forward] != 0.0) 
      d->sig[DATA::Don].weight[Signal::Forward] += spliceBoost;
    
    if ((cESTMatch & MRightForward) && 
	d->sig[DATA::Acc].weight[Signal::Forward] != 0.0)
      d->sig[DATA::Acc].weight[Signal::Forward] += spliceBoost;

    if ((cESTMatch & MLeftReverse) && 
	d->sig[DATA::Acc].weight[Signal::Reverse])
      d->sig[DATA::Acc].weight[Signal::Reverse] += spliceBoost;

    if ((cESTMatch & MRightReverse) &&
	d->sig[DATA::Don].weight[Signal::Reverse] != 0.0)
      d->sig[DATA::Don].weight[Signal::Reverse] += spliceBoost;

    // Exon Forward
    // Si on a un Gap EST ou si l'on connait le sens du match EST
    for(int i=0; i<3; i++)
      if ((cESTMatch & Gap) || 
	  ((cESTMatch & Hit) && !(cESTMatch & HitForward)))
	d->contents[i] += estP;
    
    // Exon Reverse
    // Si on a un Gap EST ou si l'on connait le sens du match EST
    for(int i=3; i<6; i++)              
      if ((cESTMatch & Gap) ||
	  ((cESTMatch & Hit) && !(cESTMatch & HitReverse)))
	d->contents[i] += estP;
    
    // Intron Forward
    // Si on a un Hit EST ou si l'on connait le sens du match EST
    if((cESTMatch & Hit) ||
       ((cESTMatch & Gap) && !(cESTMatch & GapForward))) {
      d->contents[DATA::IntronF] += estP;
      d->contents[DATA::IntronUTRF] += estP;
    }
    
    // Intron Reverse
    // Si on a un Hit EST ou si l'on connait le sens du match EST
    if((cESTMatch & Hit) ||
       ((cESTMatch & Gap) && !(cESTMatch & GapReverse))) {
      d->contents[DATA::IntronR] += estP;
      d->contents[DATA::IntronUTRR] += estP;
    }

    // UTR Forward
    // Si on a un hit ou un gap ET qu'il est sur l'autre brin seulement
    if ((cESTMatch & (Hit | Gap)) && !(cESTMatch & (HitForward | GapForward))) {
      d->contents[DATA::UTR5F] += estP;
      d->contents[DATA::UTR3F] += estP;
    }
    
    // UTR Reverse
    // Si on a un hit un un gap ET qu'il est sur l'autre brin seulement
    if ((cESTMatch & (Hit | Gap)) && !(cESTMatch & (HitReverse | GapReverse))) {
      d->contents[DATA::UTR5R] += estP;
      d->contents[DATA::UTR3R] += estP;
    }

    // Ca merdouille sur SeqAra. A creuser
    // UTR Forward
    // Si on a un hit ou un gap ET qu'il est sur l'autre brin seulement
//     if ((cESTMatch & (Hit | Gap)) && !(cESTMatch & (HitForward | GapForward))) {
//       d->contents[DATA::UTR5F] += estP;
//       d->contents[DATA::UTR3F] += estP;
//     }
    
//     // UTR Reverse
//     // Si on a un hit un un gap ET qu'il est sur l'autre brin seulement
//     if ((cESTMatch & (Hit | Gap)) && !(cESTMatch & (HitReverse | GapReverse))) {
//       d->contents[DATA::UTR5R] += estP;
//       d->contents[DATA::UTR3R] += estP;
//     }
    
    // Intergenique: tout le temps si on a un match
    d->contents[DATA::InterG] += ((cESTMatch & (Gap|Hit)) != 0)*estP;
    

    d->EstMatch = cESTMatch;  // WARNING : EST -> on est dans intron
    index++;
  }

   // Pour que les UTR soient supportés par un EST
  if (cESTMatch == 0  &&  (int)vPos.size() != 0) {
    // Left
    for (int k=1; k<=utrM; k++) {
      iter = lower_bound(vPos.begin(), vPos.end(), pos+k);
      if(*iter == pos+k) {
	cESTMatch = vESTMatch[iter-vPos.begin()];
	// If only Margin (-> EST extremities) then penalize all utr tracks
	if((cESTMatch & Margin) && !(cESTMatch & Gap) && !(cESTMatch & Hit))
	  {
	    d->contents[DATA::UTR5F] += log(utrP);
	    d->contents[DATA::UTR5R] += log(utrP);
	    d->contents[DATA::UTR3F] += log(utrP);
	    d->contents[DATA::UTR3R] += log(utrP);
	    break;
	  }
      }
    }
    // Right
    for (int k=utrM; k>0; k--) {
      iter = lower_bound(vPos.begin(), vPos.end(), pos-k);
      if(*iter == pos-k) {
	cESTMatch = vESTMatch[iter-vPos.begin()];
	// If only Margin (-> EST extremities) then penalize all utr tracks
	if((cESTMatch & Margin) && !(cESTMatch & Gap) && !(cESTMatch & Hit))
	  {
	    d->contents[DATA::UTR5F] += log(utrP);
	    d->contents[DATA::UTR5R] += log(utrP);
	    d->contents[DATA::UTR3F] += log(utrP);
	    d->contents[DATA::UTR3R] += log(utrP);
	    break;
	  }
      }
    }
  }
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
  
  fprintf(stderr,"Reading cDNA hits............");
  fflush(stderr);

  AllEST = AllEST->ReadFromFile(ESTFile,NumEST);
  
  fprintf(stderr,"%d sequences read\n",*NumEST);
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
    double WorstSpliceF, WorstSpliceR;
    double DonF,AccF,DonR,AccR;
    DATA dTmp;
    
    ThisEST = HitTable[index];

    // Le veritable  brin est a priori indetermine
    TheStrand = HitForward | HitReverse;
    Inc = 0;
    ExonInc = 0;
    WorstSpliceF = WorstSpliceR = -NINFINITY;
        
    // Look for each match in the current Hit
    ThisBlock = ThisEST->Match;
   
    // First Step: tries to determine strand
    while (ThisBlock) {
      // si ona un gap ?
      if ((ThisBlock->Prev != NULL) && 
	  abs(ThisBlock->Prev->LEnd - ThisBlock->LStart) <= 6) {
	DonF = NINFINITY;
	DonR = NINFINITY;
	AccF = NINFINITY;
	AccR = NINFINITY;
	
	for (j = -EstM; j <= EstM; j++) {
	  k = Min(X->SeqLen,Max(0,ThisBlock->Prev->End+j+1));
	  MS->GetInfoSpAt(Type_Acc|Type_Don, X, k, &dTmp);

	  DonF = Max(DonF, dTmp.sig[DATA::Don].weight[Signal::Forward]-
		     dTmp.sig[DATA::Don].weight[Signal::ForwardNo]);
	  AccR = Max(AccR, dTmp.sig[DATA::Acc].weight[Signal::Reverse]-
		     dTmp.sig[DATA::Acc].weight[Signal::ReverseNo]);
	  
	  k = Min(X->SeqLen,Max(0,ThisBlock->Start+j));
	  if(MS->GetInfoSpAt(Type_Acc|Type_Don, X, k, &dTmp)) {
	    DonR = Max(DonR, dTmp.sig[DATA::Don].weight[Signal::Reverse]-
		       dTmp.sig[DATA::Don].weight[Signal::ReverseNo]);
	    AccF = Max(AccF, dTmp.sig[DATA::Acc].weight[Signal::Forward]-
		       dTmp.sig[DATA::Acc].weight[Signal::ForwardNo]);
	  }
	  else {
	    fprintf(stderr,"   WARNING: cDNA hits ignored."
		    " No splices sites predicted !\n");
	    exit(2);
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

    //    printf("Extreme splices: %e %e\n",WorstSpliceF,WorstSpliceR);

    // Tous les blocs ont ete traites
    if (WorstSpliceF == NINFINITY) TheStrand &= (~HitForward);
    if (WorstSpliceR == NINFINITY) TheStrand &= (~HitReverse);
    
    // next iteration on the same EST
    ThisBlock = ThisEST->Match;
    while (TheStrand && ThisBlock) {
      // Check for consistency with already read Hits
      // The Inc flag will keep the Inconsistency status
      // 1 - inconsistent with a previous EST
      // 2 - no splice site on the borders of a gap
      // 3 - an exon contains STRONG donor on both strands

      for (i = ThisBlock->Start+EstM; !Inc && i <= ThisBlock->End-EstM; i++) {

	if (((ESTMatch[i] & Hit) && !(ESTMatch[i] & TheStrand)) ||
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
	  if (((ESTMatch[i] & Gap) && !(ESTMatch[i] & (TheStrand << HitToGap))) ||
	    (Inconsistent(ESTMatch[i] | (TheStrand << HitToGap)))) {
            fprintf(stderr,"   [%s]: inconsistent gap [%d-%d]\n",
                    ThisEST->Name,ThisBlock->Prev->End+2,ThisBlock->Start);
            Inc = 1;
          }
      }

      DonF = NINFINITY;
      DonR = NINFINITY;
      
      // calcul des sites d'epissage internes a l'exon
      for (i = ThisBlock->Start+EstM+1; !Inc && i <= ThisBlock->End-EstM-1; i++) {
	MS->GetInfoSpAt(Type_Acc|Type_Don, X, i, &dTmp);
	DonF = Max(DonF, dTmp.sig[DATA::Don].weight[Signal::Forward]-
		   dTmp.sig[DATA::Don].weight[Signal::ForwardNo]);
	DonR = Max(DonR, dTmp.sig[DATA::Don].weight[Signal::Reverse]-
		   dTmp.sig[DATA::Don].weight[Signal::ReverseNo]);
      }
      if (DonF > DonorThreshold) ExonInc |= 1;
      if (DonR > DonorThreshold) ExonInc |= 2;
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
      unsigned char Info;

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
	  if ((Info = (ESTMatch[i] & Hit)))
	    ESTMatch[i] |= (Info & TheStrand);
	  else
	    ESTMatch[i] |= TheStrand;
	}

	// Do we want to mark extreme region or alignements as margins ?
#ifdef EXTREMEMARGIN	
	if (ThisBlock->Prev == NULL)
	  for (i = ThisBlock->Start; i<LBoundary; i++)
	    ESTMatch[i] |= TheStrand << HitToMLeft;

	if (ThisBlock->Next == NULL)
	  for (i = RBoundary+1; i <= ThisBlock->End; i++)
	    ESTMatch[i] |= TheStrand << HitToMRight;
#endif
	
	if ((ThisBlock->Prev != NULL) && 
	    abs(ThisBlock->Prev->LEnd-ThisBlock->LStart) <= 6) {
	  
	  for (i = Max(0,ThisBlock->Prev->End+1-EstM); 
	       i < Min(X->SeqLen, ThisBlock->Prev->End+1+EstM); i++) 
	    if ((Info = (ESTMatch[i] & MLeft)))
	      ESTMatch[i] |= (Info & (TheStrand << HitToMLeft));
	    else
	      ESTMatch[i] |= TheStrand << HitToMLeft;
	  
	  for (i = ThisBlock->Prev->End+1; i <= ThisBlock->Start-1; i++)
	    if ((Info = (ESTMatch[i] & Gap)))
	      ESTMatch[i] |= (Info & (TheStrand << HitToGap));
	    else
	      ESTMatch[i] |= TheStrand << HitToGap;

	  for (i =  Max(0,ThisBlock->Start-1-EstM);
	       i <  Min(X->SeqLen, ThisBlock->Start-1+EstM); i++)
	    if ((Info = (ESTMatch[i] & MRight)))
	      ESTMatch[i] |= (Info & (TheStrand << HitToMRight));
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

// ----------------------------
//  Plot Sensor information
// ----------------------------
void SensorEst :: Plot(DNASeq *TheSeq)
{
}

// ------------------
//  Post analyse
// ------------------
void SensorEst :: PostAnalyse(Prediction *pred)
{
  int TStart = 0, TEnd = 0, GStart = 0, GEnd = 0;
  int state, stateNext = 0, stateBack = 0;
  int pos,   posNext   = 0, posBack   = 0;
  int NumGene = 0;  

  int pprocess = PAR.getI("Est.PostProcess",GetNumber());
  if (!pprocess) return;

  if (NumEST == 0) return;

  qsort((void *)HitTable,NumEST,sizeof(void *),HitsCompareLex);
  
  // Reset static in EST Support or FEA Support
  if (pprocess == 1)  ESTSupport(NULL,100,0,100,0,NULL,0);
  else                FEASupport(NULL,100,0,100,0,NULL,0,0);

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
    
    // Demarrage exon extreme. Noter pour verif EST
    if ((state     == UTR5F) || (state     == UTR3R)) GStart = pos;
    if ((stateNext == UTR3F) || (stateNext == UTR5R)) GEnd   = pos-1;
    if ((state     == UTR5F) || (state     == UTR3R)) TStart = posBack;
    if ((stateNext == UTR3F) || (stateNext == UTR5R) || i==0) TEnd = posNext-1;
        
    // Si fin de gene ou fin seq gene en cours
    if ((state <= TermR3 && stateNext >= InterGen) ||
	(i==0 && (state <= IntronR3 || state == UTR3R || state == UTR5F))) {
      NumGene++;
      if (pprocess == 1)
	// Analyse des ESTs par rapport à la prédiction
	ESTSupport(pred,TStart,TEnd,GStart,GEnd,HitTable,NumEST);
      else
	// Analyse de la prédiction (des features) par rapport aux EST
	FEASupport(pred,TStart,TEnd,GStart,GEnd,HitTable,NumEST,NumGene);
      GStart = TStart = GEnd = TEnd = -1;
    }
  }
}

// -------------------------------------------------------------------------
//  Post analyse (Analyse des ESTs par rapport à la prédiction)
//    Tdebut/fin = debut/fin de transcript
//    debut/fin  = debut/fin de traduit
// -------------------------------------------------------------------------
void SensorEst :: ESTSupport(Prediction *pred, int Tdebut, int Tfin,
			     int debut, int fin,  Hits **HitTable, int Size)
{
  static int EstIndex;
  int supported = 0;
  int CDSsupported = 0;
  unsigned char *Sup;
  Block *ThisBlock;
  int ConsistentEST,i;
  int from = 0, to = 0, ESTEnd = 0;
  int state;

  if (pred == NULL) {
    EstIndex = 0;
    return;
  }
  
  Sup = new unsigned char[Tfin-Tdebut+1]; 
  
  for (i=0; i <= Tfin-Tdebut; i++)
    Sup[i]=0;
  
  // si la fin du codant n'a pas ete rencontree
  if (fin == -1) fin = Tfin;
  if ((debut == -1) || (debut > Tfin)) debut = Tfin+1;
  
  //si l'iteration precendete a atteint l'extremite 
  if (EstIndex >= Size) EstIndex = Max(0,Size-1);

  // on rembobine....
  while ((EstIndex > 0) && (HitTable[EstIndex]->End > Tdebut))
    EstIndex--;

  if (EstIndex >= 0  &&  HitTable[EstIndex]->End < Tdebut) EstIndex++;
  
  while (EstIndex >=0  &&  EstIndex < Size) {
    // le dernier transcrit exploitable est passe
    if (HitTable[EstIndex]->Start > Tfin) break;
    
    ConsistentEST = 1;
    ThisBlock = HitTable[EstIndex]->Match;
    
    while (ThisBlock && ConsistentEST) {
      // si on a un gap
      if (ThisBlock->Prev != NULL) {
	// intersection
	from = Max(Tdebut,ThisBlock->Prev->End+1);
	to = Min(Tfin,ThisBlock->Start-1);
	
	for (i = from; i <= to; i++) {
	  state = pred->getStateForPos(i+1);
	  if ((state < IntronF1) || state == InterGen)
	    ConsistentEST = 0;
	}
      }      
      from = Max(Tdebut,ThisBlock->Start);
      to = Min(Tfin,ThisBlock->End);
      ESTEnd = ThisBlock->End;
      
      for (i = from; i <= to; i++) {
	state = pred->getStateForPos(i+1);
	if ((state > TermR3) && (state <= InterGen))
	  ConsistentEST = 0;
      }
      ThisBlock = ThisBlock->Next;
    }
    printf("cDNA  %-12s %7d %7d     %4d     %2d introns    ",
	   HitTable[EstIndex]->Name,
	   HitTable[EstIndex]->Start+1,HitTable[EstIndex]->End+1,
	   HitTable[EstIndex]->Length,HitTable[EstIndex]->NGaps);
    
    if (HitTable[EstIndex]->Rejected) printf("Filtered ");
    else if (!ConsistentEST) printf("Inconsistent");
    else {
      if ((HitTable[EstIndex]->Start) <= Tdebut && (ESTEnd >= Tfin))
	printf("Full Transcript Support");
      else if ((HitTable[EstIndex]->Start) <= debut && (ESTEnd >= fin))
	printf("Full Coding Support");
      else printf("Support");
       
      ThisBlock = HitTable[EstIndex]->Match;
      while (ThisBlock) {
	if (ThisBlock->Prev != NULL) {
	  // intersection
	  from = Max(Tdebut,ThisBlock->Prev->End+1);
	  to = Min(Tfin,ThisBlock->Start-1);
	  
	  for (i= from; i <= to; i++) 
	    if (!Sup[i-Tdebut]) {
	      Sup[i-Tdebut] = TRUE;
	      supported++;
	      if ((i >=debut) && (i <=fin)) CDSsupported++;
	    }
	}
	from = Max(Tdebut,ThisBlock->Start);
	to = Min(Tfin,ThisBlock->End);
	
	for (i= from; i <= to; i++) 
	  if (!Sup[i-Tdebut]) {
	    Sup[i-Tdebut] = TRUE;
	    supported++;
	    if ((i >=debut) && (i <=fin)) CDSsupported++;
	  }
	ThisBlock = ThisBlock->Next;
      }
    }
    printf("\n");
    EstIndex++;
  }
  if (fin >= debut)
    printf("      CDS          %7d %7d    %5d     supported on %d bases\n",
	   debut+1,fin+1,fin-debut+1,CDSsupported);
  printf("      Gene         %7d %7d    %5d     supported on %d bases\n",
	 Tdebut+1,Tfin+1,Tfin-Tdebut+1,supported);
  delete Sup;
  return;
}

// -------------------------------------------------------------------------
//  Post analyse (Analyse de la prédiction (features) par rapport aux EST)
//    Tdebut/fin = debut/fin de transcript
//    debut/fin  = debut/fin de traduit
// -------------------------------------------------------------------------
void SensorEst :: FEASupport(Prediction *pred, int Tdebut, int Tfin, int debut,
			     int fin, Hits **HitTable, int Size, int NumGene)
{
  static int EstIndex;
  Block *ThisBlock;
  int ConsistentEST, i, j;
  int from = 0, to = 0, ESTEnd = 0;
  int state;
  std::vector <int> vSupEstI;  // index des transcrits supportant la pred
  
  if (pred == NULL) {
    EstIndex = 0;
    return;
  }
  
  /***********************************************************************/
  /** Objectif : obtenir un vecteur contenant les index des transcripts **/
  /**            qui supportent la prédiction                           **/
  /***********************************************************************/
  // si la fin du codant n'a pas ete rencontree
  if (fin == -1) fin = Tfin;
  if ((debut == -1) || (debut > Tfin)) debut = Tfin+1;
  
  //si l'iteration precendete a atteint l'extremite 
  if (EstIndex >= Size) EstIndex = Max(0,Size-1);

  // on rembobine....
  while ((EstIndex > 0) && (HitTable[EstIndex]->End > Tdebut))
    EstIndex--;

  if (EstIndex >= 0  &&  HitTable[EstIndex]->End < Tdebut) EstIndex++;
  
  while (EstIndex >=0  &&  EstIndex < Size) {
     // le dernier transcrit exploitable est passe
    if (HitTable[EstIndex]->Start > Tfin) break;
    
    ConsistentEST = 1;
    ThisBlock = HitTable[EstIndex]->Match;
    
    while (ThisBlock && ConsistentEST) {
      // si on a un gap
      if (ThisBlock->Prev != NULL) {
	// intersection
	from = Max(Tdebut,ThisBlock->Prev->End+1);
	to = Min(Tfin,ThisBlock->Start-1);
	
	for (i = from; i <= to; i++) {
	  state = pred->getStateForPos(i+1);
	  if ((state < IntronF1) || state ==  InterGen )
	    ConsistentEST = 0;
	}
      }      
      from = Max(Tdebut,ThisBlock->Start);
      to = Min(Tfin,ThisBlock->End);
      ESTEnd = ThisBlock->End;
      
      for (i = from; i <= to; i++) {
	state = pred->getStateForPos(i+1);
	if ((state > TermR3) && (state <= InterGen))
	  ConsistentEST = 0;
      }
      ThisBlock = ThisBlock->Next;
    }
    
    // Si EST non filtrée et "consistente"
    if (!HitTable[EstIndex]->Rejected  &&  ConsistentEST) {
      vSupEstI.push_back( EstIndex );
    }
    EstIndex++;
  }

  /***********************************************************************/
  /* Objectif : Analyser chaque feature predite -> supportée ?           */
  /***********************************************************************/
  int  start, end, len;
  char fea[5];
  char strand = '+';
  if (PhaseAdapt(pred->getStateForPos(debut+1)) < 0) strand = '-';
  int NumFEA = 0;
  if(strand == '-') NumFEA = pred->nbExon(NumGene) + 1;
    
  for(i=pred->size()-1; i!=-1; i--) {
    state = pred->getState(i);
    start = pred->getPos(i+1) + 1;
    end   = pred->getPos(i);
    if(end >= Tdebut-1  &&  end <= Tfin+1)
      {
	len      = 0;
	int numF = -1;

	if(state <= TermR3) {
	  (strand == '-') ? NumFEA-- : NumFEA++;
	  strcpy(fea, "Exon");
	  numF = NumFEA;
	}
	else if(state == UTR5F  ||  state == UTR5R) {
	  strcpy(fea, "UTR5");
	  numF = 0;
	}
	else if(state == UTR3F  ||  state == UTR3R) {
	  strcpy(fea, "UTR3");
	  numF = 0;
	}

	if(numF != -1) {
	  if((int)vSupEstI.size() > 0)
	    // Longueur totale supportée par les transcrits
	    len = LenSup(HitTable, vSupEstI, -1, start, end);
	  
	  if(len > 0) {
	    printf("%d.%d\tEuGene_cDNA\t%s\t%d\t%d\t%d\t%c\t.\t",
		   NumGene, numF, fea, start, end, len, strand);
	    printf("%.1f\t", (float)len/(end-start+1)*100);
	    
	    for(j=0; j<(int)vSupEstI.size(); j++) {
	      // Longueur supportée par le transcrit
	      len =  LenSup(HitTable, vSupEstI, j, start, end);
	      if(len > 0)
		printf("%s(%.1f) ", HitTable[vSupEstI[j]]->Name,
		       (float)len/(end-start+1)*100);
	    }
	    printf("\n");
	  }
	}
      }	
  }
  
  /***********************************************************************/
  /* Objectif : Analyser la CDS et le gene predit -> supportée ?         */
  /***********************************************************************/
  for(int i=0; i<2; i++) {
    if(i==0) { start = debut;  end = fin;  strcpy(fea, "CDS");  }
    else     { start = Tdebut; end = Tfin; strcpy(fea, "Gene"); }

    if (end >= start) {
      len = 0;
      if((int)vSupEstI.size() > 0)
	// Longueur totale supportée par les transcrits
	len = LenSup(HitTable, vSupEstI, -1, start, end);
      
      if(len > 0) {
	printf("%d\tEuGene_cDNA\t%s\t%d\t%d\t%d\t%c\t.\t",
	       NumGene, fea, start+1, end+1, len, strand);
	printf("%.1f\t", (float)len/(end-start+1)*100);
	
	for(j=0; j<(int)vSupEstI.size(); j++) {
	  // Longueur supportée par le transcrit
	  len =  LenSup(HitTable, vSupEstI, j, start, end);
	  if(len > 0)
	    printf("%s(%.1f) ", HitTable[vSupEstI[j]]->Name,
		   (float)len/(end-start+1)*100);
	}
	printf("\n");
      }
    }
  }
  vSupEstI.clear();
  return;
}

// -------------------------------------------------------------------------
//  Lenght supported by EST.
//    If index = -1 -> Lenght supported by ALL EST.
//    Else          -> Lenght supported by ONE EST.
// -------------------------------------------------------------------------
int SensorEst :: LenSup(Hits **HitTable, std::vector<int> vSupEstI,
			int index, int beg, int end)
{
  int supported = 0;
  unsigned char *Sup;
  Block *ThisBlock;
  int from = 0, to = 0;
  int i;
  int j;

  Sup = new unsigned char[end-beg+1]; 
  for (i=0; i<=end-beg; i++)
    Sup[i]=0;
  
  index == -1  ?  i=0 : i=index;
  for(; i<(int)vSupEstI.size() || i==index; i++) {
    ThisBlock = HitTable[vSupEstI[i]]->Match;
    while (ThisBlock) {
      if (ThisBlock->Prev != NULL) {
 	// intersection
 	from = Max(beg, ThisBlock->Prev->End+1);
 	to   = Min(end, ThisBlock->Start-1);

 	for (j=from; j<=to; j++)
 	  if (!Sup[j-beg]) {
 	    Sup[j-beg] = TRUE;
 	    supported++;
 	  }
      }
      from = Max(beg, ThisBlock->Start);
      to   = Min(end, ThisBlock->End);
      
      for (j=from; j<=to; j++)
	if (!Sup[j-beg]) {
	  Sup[j-beg] = TRUE;
	  supported++;
	}
      ThisBlock = ThisBlock->Next;
    }
    if(index != -1) break;
  }
  delete Sup;
  return supported;
}
