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
    printf("   pos nt  EF1   EF2   EF3   ER1   ER2   ER3    IF    IR    IG   U5F   U5R   U3F   U3R FW: tSta tSto  Sta  Sto  Acc  Don  Ins  Del REV: tSta tSto  Sta  Sto  Acc  Don  Ins  Del\n");
    for(int i=0; i<Data_Len ; i++) {
      MS.GetInfoAt   (X, i, &Data);      MS.PrintDataAt (X, i, &Data);
    }
    printf("   pos nt  EF1   EF2   EF3   ER1   ER2   ER3    IF    IR    IG   U5F   U5R   U3F   U3R FW: tSta tSto  Sta  Sto  Acc  Don  Ins  Del REV: tSta tSto  Sta  Sto  Acc  Don  Ins  Del\n");
  }
  
  else if ((printopt0 == 'l') || (printopt0 == 'h') ||
	   (printopt0 == 'g') || (printopt0 == 'a')) {
    int nbGene  = 1;
    int nbExon  = 0;
    int cons = 0, incons = 0;
    int forward,init,term,Lend,Rend,Phase;
    int Don,Acc;
    char seqn[6] = "";
    char *position;
    int stateBack = 0, state, stateNext = 0;
    int posBack   = 0, pos;
        
    fprintf(stderr,"\n");
      
    if (printopt0 == 'h') {
      printf("<HTML><TITLE>EuGene</TITLE><BODY><CENTER><H1>EuGene prediction</H1></CENTER>\n");
      printf("<center><listing>\n");
      printf("\n\t      Type    S       Lend    Rend   Length  Phase   Frame      Ac      Do     Pr.\n");
    }
    
    else if (printopt0 == 'l')
      fprintf(stderr,"    Seq         Type    S       Lend    Rend   Length  Phase   Frame      Ac      Do     Pr.\n");
    
    else if (printopt0 == 'a')
      fprintf(stderr,"Seq   Type    S       Lend    Rend   Length  Phase   Frame      Ac      Do     Pr.\n");
    
    if (printopt0 == 'g' && sequence == optind)
      printf("name\tsource\tfeature\tstart\tend\tscore\tstrand\tframe\n");
    
    if(printopt0 != 'g')
      if(sequence != optind)
	printf("\n");
      else fprintf(stderr,"\n");
    
    //position = strstr(argv[sequence],"/seq");
    position = BaseName(argv[sequence]);
    if (position  == NULL)
      strcpy(seqn,"     ");
    else {
      if (char * suffix = rindex(position,'.')) *suffix = 0; // on enleve l'extension (.fasta)
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
      if(i != 0)
	stateNext = pred->getState(i-1);
      
      if(pos == 0)
	continue;
 
      if (estopt)
       	CheckConsistency(posBack, pos, state, &cons, &incons, X);
      
      // An exon is finishing
      if (state <= ExonR3) {
	// strand ?
	forward = (state < 3);
	if(forward) nbExon++;
	else        nbExon--;
	if(!forward && (i == pred->size()-1 ||
			(i == pred->size()-2 && stateBack < InterGen5)))
	  nbExon = pred->nbExon(1);
	
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
	
	if(printopt0 == 'g' || printopt0 == 'a')
	  printf("%s",seqn);
	else
	  printf("%s.%d.%d.%d",seqn,sequence-optind+1,nbGene,nbExon);

	if (printopt0 == 'g') printf("\tEuGene\t");
	else if (printopt0 == 'a') printf(" ");
	else printf("\t");
	
	if (init && term) {
	  printf("Sngl");
	  nbExon = 0;
	}
	else if (init) {
	  printf("Init");
	  if(!forward) nbExon = 0;
	}
	else if (term) {
	  printf("Term");
	  if(forward)  nbExon = 0;
	}
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
	if(printopt0 == 'g' || printopt0 == 'a')
	  printf("%s",seqn);
	else
	  printf("%s.%d.%d.%d",seqn,sequence-optind+1,nbGene,nbExon);
	
	if (printopt0 == 'g') printf("\tEuGene\t");
	else if (printopt0 == 'a') printf(" ");
	else printf("\t");
	
	switch (state) {
	case 13: // UTR5' F
	  if (printopt0 == 'g')
	    printf("Utr5\t%d\t%d\t0\t+\t.\n",
		   offset+posBack+1, offset+pos);
	  else printf("Utr5    +");
	  break;
	  
	case 14: // UTR 3' F
	  nbGene++;
	  if (printopt0 == 'g')
	    printf("Utr3\t%d\t%d\t0\t+\t.\n",
		   offset+posBack+1, offset+pos);
	  else printf("Utr3    +");
	  break;
	  
	case 15: // UTR5' R
	  nbGene++;
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
	if(stateNext >= ExonR1 && stateNext <= ExonR3)
	  nbExon = pred->nbExon(nbGene) + 1;
      }
      if(pos == Data_Len)
	break;
    }
    
    if (printopt0 == 'h')   {
      //position = BaseName(argv[sequence]);
      position = argv[sequence];
      strcat(position,".fasta");
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

    // y a t'l de l'info
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
