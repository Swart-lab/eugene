if (printopt == 'd')
{
  for (j = 0 ; j < 8 ; j++)
    NScore [j] = 1.0;
     
  window = ((window/2)*2)+1;
     
  for (i = 0; i < window/2; i++) {
    Fill_Score(TheSeq,IMMatrix,i,BaseScore);
    for (j = 0 ; j < 8 ; j++)
      NScore [j] += log(BaseScore[j]);
  }

  printf("#  pos FR Pr  Ph 1 Ph 2 Ph 3 Ph-1 Ph-2 Ph-3 IntF IntR SF  SR  sF  sF  sR  sR  DF   DR   AF   AR\n");
  
  for  (i = 0 ; i < Data_Len ; i++)
    {
      if (i-window/2 > 0) {
        Fill_Score(TheSeq,IMMatrix,i-1-window/2,BaseScore);
	for (j = 0 ; j < 8 ; j++)
	  NScore [j] -= log(BaseScore[j]);
      }

      if (i+window/2 < Data_Len) {
        Fill_Score(TheSeq,IMMatrix,i+window/2,BaseScore);
	for (j = 0 ; j < 8 ; j++)
	  NScore [j] += log(BaseScore[j]);
      }

      for (j = 0 ; j < 8 ; j++) Score[j] = NScore[j];
	 
      AmplifyScore(Score,normopt);
      
      printf("%6d %c%c ", offset+1+i, (*TheSeq)[i],(*TheSeq)(i));
      
      PrintPhase(Choice[i+1]);

      for (j = 0 ; j < 8 ; j++)
	printf(" %.2f", Score[j]);
      
      if (Stop[0][i+1])
	printf(" %2d ", ((i+1)%3)+1);
      else
	printf("  - ");
      
      if (Stop[1][i])
	printf(" %2d ", -((Data_Len-i)%3)-1);
      else
	printf("  - ");
      
      if (Start[0][i] > 0.0)
	printf(" %2d %4.1f",1+(i%3), Start[0][i]);
      else
	printf("  -   - ");
      
      if (Start[1][i+1] > 0.0)
	printf(" %2d %4.1f",-((Data_Len-i-1)%3)-1,Start[1][i+1]);
      else
	printf("  -   - ");
      
      if (Don[0][i] > 0.0)
	printf(" %4.1f",log(Don[0][i]));
      else
	printf("  -  ");
      
      if (Don[1][i+1] > 0.0)
	printf(" %4.1f",log(Don[1][i+1]));
      else
	printf("  -  ");
      
      if (Acc[0][i] > 0.0)
	printf(" %4.1f",log(Acc[0][i]));
      else
	printf("  -  ");
      
      if (Acc[1][i+1] > 0.0)
	printf(" %4.1f",log(Acc[1][i+1]));
      else
	printf("  -  ");
      
      printf("\n");
    }
}
else if ((printopt == 'l') || (printopt == 'h') || (printopt == 'g'))
{
  int Starts[18];
  int cons =0,incons = 0;
  int TStart = 0, GStart = 0, GEnd = 0;
  int forward,init,term,Lend,Rend,Phase;
  int Don,Acc;
  char seqn[6];
  char *pos;

  if (estopt) {    
    qsort((void *)HitTable,NumEST,sizeof(void *),HitsCompareLex);
    // reset static in EST Support
    if (estanal) ESTSupport(NULL,100,0,100,0,NULL,0);
  }

  if (printopt == 'h')
    {
      printf("<HTML><TITLE>EuGene</TITLE><BODY><CENTER><H1>EuGene prediction</H1></CENTER>\n");
      printf("<center><listing>\n");
      printf("\n      Type    S       Lend    Rend   Length  Phase   Frame      Ac      Do   Pr.\n\n");
    }

  if (printopt == 'g')
    printf("name\tsource\tfeature\tstart\tend\tscore\tstrand\tframe\n");

  // Starting condition: 0  = started,  -1 = nothing started yet
  for (j = 0; j<18; j++)
    Starts[j] = ((Choice[0] == j) ? 0 : -1);

  fprintf(stderr,"\nSeq        Type    S       Lend    Rend   Length  Phase   Frame      Ac      Do   Pr.\n\n");

  //  pos = strstr(argv[sequence],"/seq");
  pos = BaseName(argv[sequence]);
  if (pos  == NULL)
    strcpy(seqn,"          ");
  else {
    *rindex(pos,'.') = 0; // on enleve l'extension (.fasta)
    strncpy(seqn,pos,10);
    if(strlen(seqn) < 10 && printopt != 'g')
      for(i=strlen(seqn); i<10; i++)
	seqn[i] = ' ';
    seqn[strlen(seqn)] = '\0';
  }
  
  // Kludge = a non existing choice to force exon termination
  if ((Choice[Data_Len] != InterGen5) && (Choice[Data_Len] != InterGen3))
    Choice[Data_Len+1] = InterGen5;
  //  Choice[0] = 120;
  
  for (i=0; i <= Data_Len; i++) {
    if (Choice[i+1] != Choice[i]) {
      // something happens

      if (estopt)
	CheckConsistency(Starts[Choice[i]],i,Choice[i],ESTMatch,&cons,&incons);
	
      // demarrage exon extreme. Noter pour verif EST
      if ((Choice[i+1] == UTR5F) || (Choice[i+1] == UTR3R)) TStart = i;
      if ((Choice[i] == UTR5F) || (Choice[i] == UTR3R)) GStart = i;
      if ((Choice[i+1] == UTR3F) || (Choice[i+1] == UTR5R)) GEnd = i-1;

      // An exon is finishing
      if (Choice[i] <= ExonR3) {
	// strand ?
	forward = (Choice[i] < 3);
	
	// first or last exon ?
	init = ((forward  && Choice[Starts[Choice[i]]] >= InterGen5) || 
		(!forward && Choice[i+1] >= InterGen5));

	term = ((!forward  && Choice[Starts[Choice[i]]] >=InterGen5) || 
		(forward && Choice[i+1] >= InterGen5));
	
	  Lend = offset+Starts[Choice[i]]+1;
	  Rend = offset+i;

	if (forward) {
	  Don = Lend-1;
	  Acc = Rend+1;
	}
	else {
	  Acc = Lend-1;
	  Don = Rend+1;
	}

	printf("%s",seqn);
	
	if (printopt == 'g') printf("\tEuGene\t");
       	else printf(" ");
	
	if (init && term) printf("Sngl");
	else if (init) printf("Init");
	else if (term) printf("Term");
	else printf ("Intr");

	if (printopt == 'g')
	  printf("\t%d\t%d\t0\t%c\t%d\n",Lend,Rend,((forward) ? '+' : '-'),PhaseAdapt(Choice[i]));
       	else {	      
	  printf("    %c    %7d %7d",((forward) ? '+' : '-'),Lend,Rend);
	  printf("     %4d  ", Rend-Lend+1);
	
	  if (init)
	    printf("   %+2d", ((forward) ? 1: -1));
	  else {
	    Phase = ((forward) ?
		     PhaseAdapt(Choice[Starts[Choice[i]]]-6) :
		     -PhaseAdapt(Choice[i+1]-9));
	    
	    if (abs(Phase) <= 3)
	      printf("   %+2d",Phase);
	    else printf(" Unk.");
	  }
	  printf("      %+2d",PhaseAdapt(Choice[i]));
	  printf(" %7d %7d ", Don,Acc);
	  printf("  %3.0f.%-3.0f\n",100.0*(double)cons/(Rend-Lend+1),
		 100.0*(double)incons/(Rend-Lend+1));
	  Starts[Choice[i]] = -1;
	}
      }
      else if ((Choice[i] >= UTR5F) && (Choice[i] <= UTR3R)) {

	printf("%s",seqn);
	
	if (printopt == 'g') printf("\tEuGene\t");
       	else printf(" ");
	
	switch (Choice[i]) {
	case 13: // UTR5' F
	  if (printopt == 'g')
	    printf("Utr5\t%d\t%d\t0\t+\t.\n", offset+Starts[Choice[i]]+1, offset+i);
	  else printf("Utr5    +");
	  break;
	  
	case 14: // UTR 3' F
	  if (printopt == 'g')
	    printf("Utr3\t%d\t%d\t0\t+\t.\n", offset+Starts[Choice[i]]+1, offset+i);
	  else printf("Utr3    +");
	  break;
	  
	case 15: // UTR5' R
	  if (printopt == 'g')
	    printf("Utr5\t%d\t%d\t0\t-\t.\n", offset+Starts[Choice[i]]+1, offset+i);
	  else printf("Utr5    -");
	  break;
	  
	case 16:// UTR 3' R
	  if (printopt == 'g')
	    printf("Utr3\t%d\t%d\t0\t-\t.\n", offset+Starts[Choice[i]]+1, offset+i);
	  else printf("Utr3    -");
	  break;
	}
	
	if(printopt != 'g') {
	  printf("    %7d %7d", offset+Starts[Choice[i]]+1,offset+i);
	  printf("     %4d  ", i-Starts[Choice[i]]);
	  printf("   NA      NA      NA      NA ");
	  printf("  %3.0f.%-3.0f\n",100.0*(double)cons/(i-Starts[Choice[i]]),
		 100.0*(double)incons/(i-Starts[Choice[i]]));
	}
	Starts[Choice[i]] = -1;
      }

      if ((Choice[i+1] == InterGen5) || (Choice[i+1] == InterGen3))
	if (estopt && estanal) {
	  ESTSupport(Choice,TStart,i-1,GStart,GEnd,HitTable,NumEST);
	  GStart = TStart = GEnd = -1;
	}

      Starts[Choice[i+1]] = i;
    }
  }
  printf("\n");
  if (printopt == 'h')   {
    //pos = BaseName(argv[sequence]);
    pos = argv[sequence];
  
    printf("</listing></center>\n");
    printf("<a href=%s.trace>Trace</a><br>",pos);
    printf("<a href=%s.starts>NetStart F</a> <a href=%s.startsR>NetStart R</a><br>",pos,pos);
    printf("<a href=%s.splices>NetGene2 F</a> <a href=%s.splicesR>NetGene2 R</a><br>",pos,pos);
    printf("<a href=%s.spliceP>SplicePred F</a> <a href=%s.splicePR>SplicePred R</a><br>",pos,pos);
    printf("<a href=%s.est>Sim4</a><br>",pos);
    printf("<a href=%s.blastx0.html>BlastX SP</a><br>",pos);
    printf("<a href=%s.blastx1.html>BlastX PIR</a><br>",pos);
    printf("<a href=%s.blastx2.html>BlastX TrEMBL</a><br>",pos);
    
    OutputHTMLFileNames();
    printf("</body></html>\n");
  }
}
else {
  int Starts[18];
  int Intergenic,decalage;
  
  Intergenic = 0;
  decalage = ((argc == optind+1) ? offset : offset*(sequence-optind));
  
  // Kludge = an intergenic state is forced at the end
  //  Choice[Data_Len+1] = 12;
 
  for (j = 0; j<18; j++)
    Starts[j] = ((Choice[0] == j) ? 0 : -1);

  for (i = 0; i <= Data_Len; i++)   {
    if (Choice[i+1] != Choice[i]) {
      Intergenic = (Choice[i+1] >= 12);
      
      for (j = 0; j < 6; j++)  {
	
	if (IsPhaseOn(Choice[i+1],j) != IsPhaseOn(Choice[i],j)) {
	  
	  if (Starts[j] != -1)    {
	    if (Starts[j] == 0) 
	      printf("%d %d ",decalage+1, decalage+i);
	    else 
	      printf("%d %d ", decalage+Starts[j]+1, decalage+i);

	    Starts[j] = -1;
	  }
	  else
	    Starts[j] = i;
	}
      }
      if (Intergenic) printf("\n");
    }
  }
}
