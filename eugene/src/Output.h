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
  int Starts[18];
  int cons =0,incons = 0;
  int TStart = 0, TEnd;
  int forward,init,term,Lend,Rend,Phase;
  int Don,Acc;
  char seqn[6];
  char *pos;

  if (estopt) {    
    qsort((void *)HitTable,NumEST,sizeof(void *),HitsCompareLex);
    // reset static in EST Support
    if (estanal) ESTSupport(0,100,0,NULL,0);
  }


  if (printopt == 'h')
    {
      printf("<HTML><TITLE>EuGene</TITLE><BODY><CENTER><H1>EuGene prediction</H1></CENTER>\n");
      printf("<center><listing>\n");
      printf("\n      Type    S       Lend    Rend   Length  Phase   Frame      Ac      Do   Pr.\n\n");
    }
  
  // Starting condition: 0  = started,  -1 = nothing started yet
  for (j = 0; j<18; j++)
    Starts[j] = ((Choice[0] == j) ? 0 : -1);

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
  if ((Choice[Data_Len] != InterGen5) && (Choice[Data_Len] != InterGen3))
    Choice[Data_Len+1] = InterGen5;
  //  Choice[0] = 120;
  
  for (i=0; i <= Data_Len; i++) {
    if (Choice[i+1] != Choice[i]) {
      // something happens

      CheckConsistency(Starts[Choice[i]],i,Choice[i],ESTMatch,&cons,&incons);
	
      // demarrage exon extreme. Noter pour verif EST
      if ((Choice[i] == UTR5F) || (Choice[i] == UTR3R)) TStart = i;
      if ((Choice[i+1] == UTR3F) || (Choice[i+1] == UTR5R)) TEnd = i-1;

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
	printf(" %7d %7d ", Don,Acc);
	printf("  %3.0f.%-3.0f\n",100.0*(double)cons/(Rend-Lend+1),
	       100.0*(double)incons/(Rend-Lend+1));
	Starts[Choice[i]] = -1;
      } 
      else if ((Choice[i] >= UTR5F) && (Choice[i] <= UTR3R)) {

	printf("%s ",seqn);

	switch (Choice[i]) {
	case 13: // UTR5' F
	  printf("Utr5    +");
	  break;

	case 14: // UTR 3' F
	  printf("Utr3    +");
	  break;

	case 15: // UTR5' R
	  printf("Utr5    -");
	  break;

	case 16:// UTR 3' R
	  printf("Utr3    -");
	  break;
	}

	printf("    %7d %7d", offset+Starts[Choice[i]]+1,offset+i);
	printf("     %4d  ", i-Starts[Choice[i]]);
	printf("   NA      NA      NA      NA ");
	printf("  %3.0f.%-3.0f\n",100.0*(double)cons/(i-Starts[Choice[i]]),
	       100.0*(double)incons/(i-Starts[Choice[i]]));
	Starts[Choice[i]] = -1;
      }

      if ((Choice[i+1] == InterGen5) || (Choice[i+1] == InterGen3))
	if (estopt && estanal) ESTSupport(Choice,TStart,TEnd,HitTable,NumEST);

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
