#include "Output.h"

extern Parameters   PAR;
extern MasterSensor MS;

void Output (DNASeq *X, Prediction *pred, int sequence, int argc, char * argv[])
{
  int  i;
  DATA Data;
  int  Data_Len  = X->SeqLen;
  char printopt0 = PAR.getC("Output.format")[0];
  int  offset    = PAR.getI("Output.offset");
  int  estopt    = PAR.getI("Sensor.Est.use");
  
  if (printopt0 == 'd') {
    MS.ResetIterator();
    for(int i=0; i<Data_Len ; i++) {
      MS.GetInfoAt   (X, i, &Data);
      MS.PrintDataAt (X, i, &Data);
    }
  }

  else if ((printopt0 == 'l') || (printopt0 == 'h') || (printopt0 == 'g')) {
    int cons = 0, incons = 0;
    int TStart = 0, GStart = 0, GEnd = 0;
    int forward,init,term,Lend,Rend,Phase;
    int Don,Acc;
    char seqn[6] = "";
    char *position;
    int stateBack = 0, state, stateNext = 0;
    int posBack   = 0, pos,   posNext   = 0;

    //if (estopt) {    
    //qsort((void *)HitTable,NumEST,sizeof(void *),HitsCompareLex);
    // reset static in EST Support
    //if (PAR.estanal) ESTSupport(NULL,100,0,100,0,NULL,0);
    //}

    if (printopt0 == 'h') {
      printf("<HTML><TITLE>EuGene</TITLE><BODY><CENTER><H1>EuGene prediction</H1></CENTER>\n");
      printf("<center><listing>\n");
      printf("\n    Type    S       Lend    Rend   Length  Phase   Frame      Ac      Do     Pr.\n\n");
    }

    fprintf(stderr,"\nSeq   Type    S       Lend    Rend   Length  Phase   Frame      Ac      Do     Pr.\n\n");
    
    if (printopt0 == 'g' && sequence == optind)
      printf("name\tsource\tfeature\tstart\tend\tscore\tstrand\tframe\n");
    
    //position = strstr(argv[sequence],"/seq");
    position = BaseName(argv[sequence]);
    if (position  == NULL)
      strcpy(seqn,"     ");
    else {
      *rindex(position,'.') = 0; // on enleve l'extension (.fasta)
      strncpy(seqn,position,5);
      if(strlen(seqn) < 5 && printopt0 != 'g')
	for(i=strlen(seqn); i<5; i++)
	  strcat(seqn, " ");
      seqn[5] = '\0';
    }
       
    for(i=pred->size()-1; i!=-1; i--) {
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
      
      if (estopt)
       	CheckConsistency(posBack, pos, state, &cons, &incons, X);
      // Demarrage exon extreme. Noter pour verif EST
      if ((stateNext == UTR5F) || (stateNext == UTR3R)) TStart = pos;
      if ((state     == UTR5F) || (state     == UTR3R)) GStart = pos;
      if ((stateNext == UTR3F) || (stateNext == UTR5R)) GEnd   = pos-1;
      
      // An exon is finishing
      if (state <= ExonR3) {
	// strand ?
	forward = (state < 3);
	
	// first or last exon ?
	init = ((forward  && stateBack >= InterGen5) ||
		(!forward && stateNext >= InterGen5));

	term = ((!forward && stateBack >=InterGen5) ||
		(forward  && stateNext >= InterGen5));
	
	Lend = offset+posBack+1;
	Rend = offset+pos;
	
	if (forward) {
	  Don = Lend-1;
	  Acc = Rend+1;
	}
	else {
	  Acc = Lend-1;
	  Don = Rend+1;
	}
	
	printf("%s",seqn);
	
	if (printopt0 == 'g') printf("\tEuGene\t");
	else printf(" ");
	
	if (init && term) printf("Sngl");
	else if (init) printf("Init");
	else if (term) printf("Term");
	else printf ("Intr");
	
	if (printopt0 == 'g')
	  printf("\t%d\t%d\t0\t%c\t%d\n",
		 Lend,Rend,((forward) ? '+' : '-'),abs(PhaseAdapt(state)));
	else {
	  printf("    %c    %7d %7d",((forward) ? '+' : '-'),Lend,Rend);
	  printf("     %4d  ", Rend-Lend+1);
	  
	  if (init)
	    printf("   %+2d", ((forward) ? 1: -1));
	  else {
	    Phase = ((forward) ?
		     PhaseAdapt(stateBack-6) :
		     -PhaseAdapt(stateNext-9));
	    
	    if (abs(Phase) <= 3)
	      printf("   %+2d",Phase);
	    else printf(" Unk.");
	  }
	  printf("      %+2d",PhaseAdapt(state));
	  printf(" %7d %7d ", Don,Acc);
	  printf("  %3.0f.%-3.0f\n",100.0*(double)cons/(Rend-Lend+1),
		 100.0*(double)incons/(Rend-Lend+1));
	}
      }
      else if ((state >= UTR5F) && (state <= UTR3R)) {
	
	printf("%s",seqn);
	
	if (printopt0 == 'g') printf("\tEuGene\t");
	else printf(" ");
	
	switch (state) {
	case 13: // UTR5' F
	  if (printopt0 == 'g')
	    printf("Utr5\t%d\t%d\t0\t+\t.\n",
		   offset+posBack+1, offset+pos);
	  else printf("Utr5    +");
	  break;
	  
	case 14: // UTR 3' F
	  if (printopt0 == 'g')
	    printf("Utr3\t%d\t%d\t0\t+\t.\n",
		   offset+posBack+1, offset+pos);
	  else printf("Utr3    +");
	  break;
	  
	case 15: // UTR5' R
	  if (printopt0 == 'g')
	    printf("Utr5\t%d\t%d\t0\t-\t.\n",
		   offset+posBack+1, offset+pos);
	  else printf("Utr5    -");
	  break;
	  
	case 16:// UTR 3' R
	  if (printopt0 == 'g')
	    printf("Utr3\t%d\t%d\t0\t-\t.\n",
		   offset+posBack+1, offset+pos);
	  else printf("Utr3    -");
	  break;
	}
	
	if(printopt0 != 'g') {
	  printf("    %7d %7d", offset+posBack+1, offset+pos);
	  printf("     %4d  ",  pos - (posBack+1) +1);
	  printf("   NA      NA      NA      NA ");
	  printf("  %3.0f.%-3.0f\n",100.0*(double)cons/((offset+pos) - (offset+posBack+1)+1),
		 100.0*(double)incons/((offset+pos) - (offset+posBack+1)+1));
	}
     	//if ((Choice[i+1] == InterGen5) || (Choice[i+1] == InterGen3))
	//if (estopt && PAR.estanal) {
	//  ESTSupport(Choice,TStart,i-1,GStart,GEnd,HitTable,NumEST);
	//  GStart = TStart = GEnd = -1;
	//}
      }
    }
    
    if (printopt0 != 'g')  printf("\n");
    
    if (printopt0 == 'h')   {
      //position = BaseName(argv[sequence]);
      position = argv[sequence];
      
      printf("</listing></center>\n");
      printf("<a href=%s.trace>Trace</a><br>",position);
      printf("<a href=%s.starts>NetStart F</a> <a href=%s.startsR>NetStart R</a><br>",position,position);
      printf("<a href=%s.splices>NetGene2 F</a> <a href=%s.splicesR>NetGene2 R</a><br>",position,position);
      printf("<a href=%s.spliceP>SplicePred F</a> <a href=%s.splicePR>SplicePred R</a><br>",position,position);
      printf("<a href=%s.est>Sim4</a><br>",position);
      printf("<a href=%s.blastx0.html>BlastX SP</a><br>",position);
      printf("<a href=%s.blastx1.html>BlastX PIR</a><br>",position);
      printf("<a href=%s.blastx2.html>BlastX TrEMBL</a><br>",position);
      
      OutputHTMLFileNames();
      printf("</body></html>\n");
    }
  }
  
  else {
    int decalage;
    int line = 0;
    int state, posBack = 0, pos;
    
    decalage = ((argc == optind+1) ? offset : offset*(sequence-optind));
    
    for(i=pred->size()-1; i!=-1; i--) {
      if(i != pred->size()-1)
	posBack   = pred->getPos(i+1);
      state = pred->getState(i);
      pos   = pred->getPos(i);
      
      if (state <= ExonR3) {
	printf("%d %d ",decalage+posBack+1, decalage+pos);
	line = 0;
      }
      else if(i != pred->size()-1 && state == InterGen5 || state == InterGen3)
	{
	  line = 1;
	  printf("\n");
	}
    }
    if(line)
      printf("\n");
    else
      printf("\n\n");
  }
}

//--------------------------------------------------
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

//--------------------------------------------------
int IsPhaseOn(char p, int pp)
{
  if (p == 6) return FALSE;
  if (p == pp) return TRUE;
  return FALSE;
}

//--------------------------------------------------
// Convertit les phases 0-6 en 1 2 3 -1 -2 -3 0
//--------------------------------------------------
int PhaseAdapt(char p)
{
  if (p >= 12) return 0;
  else if (p < 3) return (1+p);
  else if (p < 6) return (2-p);
  else if (p < 9) return (p-2);
  else return (5-p);
}

//--------------------------------------------------
void PrintPhase(char p)
{
  switch (p) {
  case 0:
    printf("E1 ");
    return;
    
  case 1:
    printf("E2 ");
    return;

  case 2:
    printf("E3 ");
    return;

  case 3:
    printf("e1 ");
    return;

  case 4:
    printf("e2 ");
    return;

  case 5:
    printf("e3 ");
    return;

  case 6:
    printf("I1 ");
    return;

  case 7:
    printf("I2 ");
    return;

  case 8:
    printf("I3 ");
    return;

  case 9:
    printf("i1 ");
    return;

  case 10:
    printf("i2 ");
    return;

  case 11:
    printf("i3 ");
    return;

  case 12:
    printf("IG ");
    return;

  case 13:
    printf("U5 ");
    return;

  case 14:
    printf("U3 ");
    return;

  case 15:
    printf("u5 ");
    return;

  case 16:
    printf("u3 ");
    return;

  case 17:
    printf("IG ");
    return;

  default:
    fprintf(stderr,"ERROR: unexpected Choice value\n");
    exit(1);
  }
  return;
}

//--------------------------------------------------
// Dump les signaux au format util user
//--------------------------------------------------
void DumpSignals(int Len, REAL** Donor, REAL **Acceptor,REAL** Stop, REAL **Start,FILE* flot)
{
  int i;
  for (i=0; i< Len; i++) {
    //    if (Stop[0][i]) fprintf(flot,"stop f %d %a\n",i+1,Stop[0][i]);
    //    if (Stop[1][i]) fprintf(flot,"stop r %d %a\n",i+1,Stop[1][i]);

    if (Start[0][i]) fprintf(flot,"start f %d %a\n",i+1,Start[0][i]);
    if (Start[1][i]) fprintf(flot,"start r %d %a\n",i+1,Start[1][i]);

    if (Donor[0][i]) fprintf(flot,"donor f %d %a\n",i+1,Donor[0][i]);
    if (Donor[1][i]) fprintf(flot,"donor r %d %a\n",i+1,Donor[1][i]);

    if (Acceptor[0][i]) fprintf(flot,"acceptor f %d %a\n",i+1,Acceptor[0][i]);
    if (Acceptor[1][i]) fprintf(flot,"acceptor r %d %a\n",i+1,Acceptor[1][i]);
  }
}

// -------------------------------------------------------------------------
// Verif coherence EST: calcul le nombre de nuc. coherents et
// incoherents avec les match est
// debut/fin/etat: debut et fin de la seq. dont l'etat est etat
// cons/incons: retour des valeurs
// WARNING : A modifier, utilise ESTMATCH_TMP (Cf struct DATA) !!!!!!!!!!!!
// -------------------------------------------------------------------------
void CheckConsistency(int debut, int fin, int etat, 
		      int * cons, int* incons, DNASeq *X)
{
  int i, con = 0, inc = 0;
  DATA dTMP;

  // les valeurs qui sont coherentes avec chaque etat
  const unsigned char Consistent[18] = {
    HitForward|MForward,    HitForward|MForward,    HitForward|MForward,
    HitReverse|MReverse,    HitReverse|MReverse,    HitReverse|MReverse,
    GapForward|MForward,    GapForward|MForward,    GapForward|MForward,
    GapReverse|MReverse,    GapReverse|MReverse,    GapReverse|MReverse,
    0,
    HitForward|MForward|GapForward, HitForward|MForward|GapForward,
    HitReverse|MReverse|GapReverse, HitReverse|MReverse|GapReverse,
    0
  };
  
  const unsigned char MaskConsistent[18] = {
    Hit|Margin,    Hit|Margin,    Hit|Margin,    
    Hit|Margin,    Hit|Margin,    Hit|Margin,    
    Gap|Margin,    Gap|Margin,    Gap|Margin,
    Gap|Margin,    Gap|Margin,    Gap|Margin,
    0,
    Margin|Hit|Gap,    Margin|Hit|Gap,
    Margin|Hit|Gap,    Margin|Hit|Gap,
    0
  };

  if (debut == -1) debut = 0;
  
  for (i = debut; i <fin; i++) {
    
    MS.GetInfoSpAt(Type_Content, X, i, &dTMP);
    
    // y a t'il de l'info
    if (dTMP.ESTMATCH_TMP) {
      // y a t'il une info incoherente avec l'etat
      if (dTMP.ESTMATCH_TMP & ~MaskConsistent[etat]) 
	inc++;
      else if (dTMP.ESTMATCH_TMP & Consistent[etat]) 
	con++;
      else 
	inc++;
    }
  }
  *cons = con;
  *incons = inc;
}

// -------------------------------------------------------------------------
// Tdebut/fin = debut/fin de transcript
// debut/fin = debut/fin de traduit
// -------------------------------------------------------------------------
void ESTSupport(char * Choice, int Tdebut, int Tfin, int debut,int fin,  Hits **HitTable, int Size)
{
  static int EstIndex;
  int supported = 0;
  int CDSsupported = 0;
  unsigned char *Sup;
  Block *ThisBlock;
  int ConsistentEST,i;
  int from = 0, to = 0, ESTEnd = 0;
  
  if (Choice == NULL) { 
    EstIndex = 0;
    return;
  }

  Sup = new unsigned char[Tfin-Tdebut+1]; 

  for (i=0; i <= Tfin-Tdebut; i++)
    Sup[i]=0;
  
  // si la fin du codant n'a pas ete rencontree
  if (fin == -1) fin = Tfin;
  if ((debut == -1) || (debut > Tfin)) debut = Tfin+1;

  while (EstIndex < Size) {
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
	 
	 for (i = from; i <=to; i++)
	   if ((Choice[i+1] < IntronF1) || Choice[i+1] == InterGen5 || Choice[i+1] == InterGen3)
	     ConsistentEST = 0;
       }

       from = Max(Tdebut,ThisBlock->Start);
       to = Min(Tfin,ThisBlock->End);
       ESTEnd = ThisBlock->End;


       for (i = from; i <= to; i++)
	 if (((Choice[i+1] > ExonR3) && (Choice[i+1] <= InterGen5)) || Choice[i+1] == InterGen3)
	   ConsistentEST = 0;

       ThisBlock = ThisBlock->Next;
     }
      

     printf("cDNA  %-12s %7d %7d     %4d     %2d introns    ",HitTable[EstIndex]->Name,
	    HitTable[EstIndex]->Start+1,ESTEnd+1,
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
