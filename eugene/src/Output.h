if (printopt == 'd')
{
  for (j = 0 ; j < 8 ; j++)
    NScore [j] = 1.0;
     
  window = ((window/2)*2)+1;
     
  for (j = 0 ; j < 8 ; j++)
    for (i = 0; i < window/2; i++)
      NScore [j] += log(BaseScore[j][i]);

  printf("#  pos FR Pr  Ph 1 Ph 2 Ph 3 Ph-1 Ph-2 Ph-3 IntF IntR SF  SR  sF  sF  sR  sR  DF   DR   AF   AR\n");
  
  for  (i = 0 ; i < Data_Len ; i++)
    {
      if (i-window/2 > 0)
	for (j = 0 ; j < 8 ; j++)
	  NScore [j] -= log(BaseScore[j][i-1-window/2]);
	 
      if (i+window/2 < Data_Len)
	for (j = 0 ; j < 8 ; j++)
	  NScore [j] += log(BaseScore[j][i+window/2]);
	 
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
else if ((printopt == 'l') || (printopt == 'h'))
{
  int Starts[6];
  int forward,init,term,Lend,Rend,Phase;
  int Don,Acc;
  char seqn[6];
  char *pos;

  if (printopt == 'h')
    {
      printf("<HTML><TITLE>EuGene</TITLE><BODY><CENTER><H1>EuGene prediction</H1></CENTER>\n");
      printf("<center><listing>\n");
      printf("\n      Type    S       Lend    Rend   Length  Phase   Frame      Ac      Do   Pr.\n\n");
    }
  
  // Starting condition: 0  = the exon is only partial
  //                     -1 = nothing started yet
  for (j = 0; j<6; j++)
    Starts[j] = (IsPhaseOn(Choice[0],j) ? 0 : -1);

  /*
  fprintf(stderr,"\n Type    S       Lend    Rend   Length  Phase   Frame      Ac      Do   Pr.\n\n");
  */
  
  fprintf(stderr,"\nSeq   Type    S       Lend    Rend   Length  Phase   Frame      Ac      Do   Pr.\n\n");

  
  pos = strstr(argv[sequence],"/seq");
  if (pos  == NULL)
    strcpy(seqn,"     ");
  else
    {
      pos++;
      strncpy(seqn,pos,5);
      seqn[5] = '\0';
    }
  
  // Kludge = a non existing choice to force exon termination
  //  Choice[Data_Len+1] = 120;
  //  Choice[0] = 120;
  
  for (i=0; i <= Data_Len; i++) {
    if (Choice[i+1] != Choice[i]) {
      // something happens
      if (Choice[i] < 6) {
	forward = (Choice[i] < 3);
	
	init = ((forward  && Choice[Starts[Choice[i]]] >= 12) || (!forward && Choice[i+1] >= 12));
	term = ((!forward  && Choice[Starts[Choice[i]]] >= 12) || (forward && Choice[i+1] >= 12));
	
	if (forward) {
	  Lend = offset+Starts[Choice[i]]+1;
	  Rend = offset+i;
	  Don = Lend-1;
	  Acc = Rend+1;
	}
	else {
	  Lend = offset+Starts[Choice[i]]+1;
	  Rend = offset+i;
	  Acc = Lend-1;
	  Don = Rend+1;
	}

	printf("%s ",seqn);
	
	if (init && term) printf("Sngl");
	else if (init) printf("Init");
	else if (term) printf("Term");
	else printf ("Intr");
	      
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
	printf(" %7d %7d   1.0", Don,Acc);
	
	printf("\n");
	Starts[Choice[i]] = -1;
      }
      
      if (Choice[i+1] < 6)
	Starts[Choice[i+1]] = i;
    }
  }
  printf("\n");
  if (printopt == 'h')   {
    pos = BaseName(argv[sequence]);

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
  int Starts[6];
  int Intergenic,decalage;
  
  Intergenic = 0;
  decalage = ((argc == optind+1) ? offset : offset*(sequence-optind));
  //  printf("\n %s \n",argv[sequence]);
  
  // Kludge = an intergenic state is forced at the end
  Choice[Data_Len+1] = 12;
 
  for (j = 0; j < 6; j++)
    Starts[j] = (IsPhaseOn(Choice[0],j) ? 0 : -1);
  
  for (i=1; i<= Data_Len; i++)   {
    if (Choice[i+1] != Choice[i]) {
      Intergenic = (Choice[i+1] >= 12);
      
      for (j = 0; j<6; j++)  {
	if (IsPhaseOn(Choice[i+1],j) != IsPhaseOn(Choice[i],j)) {
	  if (Starts[j] != -1)    {
	    if (Starts[j] == 0) {
	      if (j < 3)
		printf("%d %d ",decalage+1, decalage+i);
	      else
		printf("%d %d ",decalage+1, decalage+i);
	    }
	    else {
	      if (j < 3)
		printf("%d %d ", decalage+Starts[j]+1, decalage+i);
	      else
		printf("%d %d ", decalage+Starts[j]+1, decalage+i);
	    }
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
