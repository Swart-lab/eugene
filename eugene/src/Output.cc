#include "Output.h"
#include<iostream>

extern Parameters PAR;

std::vector <std::string>  vhtml; // (-ph)

void Output (DNASeq *X, MasterSensor* ms, Prediction *pred, int sequence, int argc, char * argv[], FILE* f)
{
  int  i;
  DATA Data;
  int  Data_Len  = X->SeqLen;
  char printopt0 = PAR.getC("Output.format")[0]; 
  int  offset    = PAR.getI("Output.offset");
  int  estopt    = PAR.getI("Sensor.Est.use");
  int  trunclen  = PAR.getI("Output.truncate");
  char nameformat[20];

  if (trunclen) sprintf(nameformat,"%%%d.%ds",trunclen,trunclen);
  else strcpy(nameformat,"%s");

  if (printopt0 == 'h') {
    if (sequence == optind) {   // -ph && first seq
      StartHTML();
      OutputHTMLFileNames(1);
      vhtml.push_back("				  </select>\n"
		      "				</td>\n"
		      "				<td width=\"100\" "
		      "align=\"center\"\n"
		      "				    bgcolor=\"#c0dbe2\">\n"
		      "				  <img src=\"WEB/Images/"
		      "next.jpg\"\n"
		      "				       onclick=\"next();\"\n"
		      "				       title=\"Next\" "
		      "align=\"middle\">\n"
		      "				  &nbsp; &nbsp; &nbsp;\n"
		      "				  <img src=\"WEB/Images/"
		      "last.jpg\"\n"
		      "				       onclick=\"last();\"\n"
		      "				       title=\"Jump to end\" "
		      "align=\"middle\">\n");
      //<!--			
      //  <input type="button" name="slidebutton"
      //     onClick="ap(this.value);"
      //     value="Start" title="AutoPlay"
      //     style="width:75;border:1 SOLID #e6e6e6;">
      //-->
      vhtml.push_back("				</td>\n"
		      "			      </tr>\n"
		      "			      <tr>\n"
		      "				<td bgcolor=\"#9bc7d0\"></td>\n"
		      "				<td align=\"center\" "
		      "bgcolor=\"white\">\n"
		      "				  <div class=\"conteneur\">\n"
		      "				    <div class=\"frame\">\n"
		      "				      <table width=\"100%%\" "
		      "cellpadding=\"1\" cellspacing=\"1\">\n"
		      "					  <tr>\n"
		      "					    <td>\n"
		      "					      <table "
		      "width=\"100%%\""
		      " cellpadding=\"2\" cellspacing=\"1\">\n");
    }
    else OutputHTMLFileNames(0);
  }
  
  if (printopt0 == 'd') {
    printf("   pos nt  EF1   EF2   EF3   ER1   ER2   ER3    IF    IR    IG   U5F   U5R   U3F   U3R   IUF   IUR FWD tSta  tSto   Sta   Sto   Acc   Don   Ins   Del REV tSta  tSto   Sta   Sto   Acc   Don   Ins   Del noF tSta  tSto   Sta   Sto   Acc   Don   Ins   Del noR tSta  tSto   Sta   Sto   Acc   Don   Ins   Del\n");
    for(int i=0; i<Data_Len ; i++) {
      ms->GetInfoAt   (X, i, &Data);      
      ms->PrintDataAt (X, i, &Data);
    }
    printf("   pos nt  EF1   EF2   EF3   ER1   ER2   ER3    IF    IR    IG   U5F   U5R   U3F   U3R   IUF   IUR FWD tSta  tSto   Sta   Sto   Acc   Don   Ins   Del REV tSta  tSto   Sta   Sto   Acc   Don   Ins   Del noF tSta  tSto   Sta   Sto   Acc   Don   Ins   Del noR tSta  tSto   Sta   Sto   Acc   Don   Ins   Del\n");
  } else 
    if ((printopt0 == 'l') || (printopt0 == 'a') || (printopt0 == 'g') ||
	(printopt0 == 'h') || (printopt0 == 'H')) { 
      int nbGene  = 1;
      int nbExon  = 0;
      int cons = 0, incons = 0;
      int forward,init,term,Lend,Rend,Phase;
      int Don,Acc;
      char *position = "";
      int stateBack = 0, state, stateNext = 0;
      int posBack   = 0, pos;
      
      int lc = 1; // EuGeneHom -pH : alternance couleur des lignes de la table
      
      if ((printopt0 != 'H') && (f==stdout))
	fprintf(stderr,"\n");
      
      
      if (printopt0 == 'h') {
	vhtml.push_back("<tr align=\"center\" class=\"fonce\">\n"
			" <td><font color=\"white\" face=\"monospace\">\n"
			"  <b>SeqName_GeneN.ExonN</b></font></td>\n"
			" <td><font color=\"white\" face=\"monospace\">\n"
			"  <b>Source</b></font></td>\n"
			" <td><font color=\"white\" face=\"monospace\">\n"
			"  <b>Feature</b></font></td>\n"
			" <td><font color=\"white\" face=\"monospace\">\n"
			"  <b>Start</b></font></td>\n"
			" <td><font color=\"white\" face=\"monospace\">\n"
			"  <b>End</b></font></td>\n"
			" <td><font color=\"white\" face=\"monospace\">\n"
			"  <b>Score</b></font></td>\n"
			" <td><font color=\"white\" face=\"monospace\">\n"
			"  <b>Strand</b></font></td>\n"
			" <td><font color=\"white\" face=\"monospace\">\n"
			"  <b>Frame</b></font></td>\n"
			"</tr>\n");
      } else 
	if (printopt0 == 'H') {
	  printf("<CENTER><LISTING><H2>EuGeneHom prediction</H2>\n");
	  printf("<TABLE width=700 CELLSPACING=2 CELLPADDING=2 class=clair>\n"
		 " <TR class=clair><TD>\n"
		 "  <TABLE CELLSPACING=2 CELLPADDING=2 border=0 width=100%%"
		 "  class=clair>\n"
		 "   <TR class=fonce>\n"
		 "    <TD width=11%% align=center><font color=white>\n"
		 "     <b>Gene<br>number</b></font></TD>\n"
		 "    <TD width=11%% align=center><font color=white>\n"
		 "     <b>Element<br>number</b></font></TD>\n"
		 "    <TD width=11%% align=center><font color=white>\n"
		 "     <b>Exons/UTR</b></font></TD>\n"
		 "    <TD width=11%% align=center><font color=white>\n"
		 "     <b>Strand</b></font></TD>\n"
		 "    <TD width=11%% align=center><font color=white>\n"
		 "     <b>Left&nbsp;end</b></font></TD>\n"
		 "    <TD width=11%% align=center><font color=white>\n"
		 "     <b>Right&nbsp;end</b></font></TD>\n"
		 "    <TD width=11%% align=center><font color=white>\n"
		 "     <b>Length</b></font></TD>\n"
		 "    <TD width=11%% align=center><font color=white>\n"
		 "     <b>Phase</b></font></TD>\n"
		 "    <TD width=11%% align=center><font color=white>\n"
		 "     <b>Frame</b></font></TD>\n"
		 "   </TR>\n");
	} else 
	  if (printopt0 == 'l') {
	    if (f==stdout)
	      fprintf(stderr,"    Seq         Type    S       Lend    Rend   Length  Phase   Frame      Ac      Do      Pr\n");
	  } else 
	    if (printopt0 == 'a') fprintf(stderr,"Seq   Type    S       Lend    Rend   Length  Phase   Frame      Ac      Do      Pr\n");
      
      if (printopt0 == 'g' && sequence == optind)
	fprintf(stderr,"name\tsource\tfeature\tstart\tend\tscore\tstrand\tframe\n");
      
      if (f==stdout) {
	if(printopt0 != 'g' && printopt0 != 'h') {
	  if (sequence != optind) printf("\n");
	  else  fprintf(stderr,"\n");
	}
      }
      
      if (printopt0 == 'H') {
	i=pred->size()-1;
	if (i == 0)
	  printf("<tr class=clair><td colspan=9 align=center><b>No exons/genes"
		 "predicted in your submitted sequence !</b></td></tr>");
      }
      if (printopt0 == 'h') {
	i=pred->size()-1;
	if (i == 0)
	  vhtml.push_back("<tr class=\"A0\">\n <td colspan=\"8\" "
			  "align=\"center\">\n  <font face=\"monospace\">\n  "
			  "<b>No exons/genes predicted in your submitted "
			  "sequence !</b></td>\n</tr>\n");
      }
      for(i=pred->size()-1; i!=-1; i--) {
	if(i != pred->size()-1) {
	  stateBack = pred->getState(i+1);
	  posBack   = pred->getPos(i+1);
	}
	state = pred->getState(i);
	pos   = pred->getPos(i);
	stateNext = (i ? pred->getState(i-1) : -1);
	
	if(pos == 0)
	  continue;
	
	if (estopt)
	  CheckConsistency(posBack, pos, state, &cons, &incons, X, ms);
	
	// An exon is finishing
	if (state <= TermR3) {
	  // strand ?
	  forward = (1 <= PhaseAdapt(state) && PhaseAdapt(state) <= 3);
	  if(forward) nbExon++;
	  else        nbExon--;
	  if(!forward && (i == pred->size()-1 ||
			  (i == pred->size()-2 && stateBack < InterGen)))
	    nbExon = pred->nbExon(1);
	  
	  // first or last exon ?
	  init = ((forward  && stateBack >= InterGen) ||
		  (!forward && stateNext >= InterGen));
	  
	  term = ((!forward && stateBack >=InterGen) ||
		  (forward  && stateNext >= InterGen));
	  
	  Lend = offset+posBack+1;
	  Rend = offset+pos;
	  
	  if (forward) {
	    Don = Lend-1;
	    Acc = Rend+1;
	  } else {
	    Acc = Lend-1;
	    Don = Rend+1;
	  }
	  
	  if(printopt0 == 'h') {
	    lc++;
	    vhtml.push_back("<tr class=\"A" + to_string(lc%2) + "\">\n"
			    " <td align=\"center\">\n"
			    "  <font face=\"monospace\">" + position + "_"
			    + to_string(nbGene) + "." + to_string(nbExon)
			    + "</font></td>\n");
	    vhtml.push_back(" <td align=\"center\">\n"
			    "  <font face=\"monospace\">");
	    vhtml.push_back(argv[0]);
	    vhtml.push_back("</font></td>\n");
	  }
	  else
	    if(printopt0 == 'g' || printopt0 == 'a')
	      fprintf(f,nameformat,X->Name);
	    else
	      if(printopt0 == 'H') {
		lc++;
		printf("   <TR class=A%d align=center>\n"
		       "    <TD>%d</TD><TD>%d</TD>",lc%2,nbGene,nbExon);
	      }
	      else {
		fprintf(f,nameformat,X->Name);
		fprintf(f,".%d.%d",nbGene,nbExon);
	      }
	  
	  if (printopt0 == 'g') 
	    printf("\tEuGene\t");
	  else 
	    if (printopt0 == 'a') 
	      printf(" ");
	    else 
	      if (printopt0 != 'H'  &&  printopt0 != 'h') fprintf(f,"\t");
	  
	  if (init && term) {
	    if (printopt0 == 'h')
	      vhtml.push_back(" <td align=\"center\">\n  <font face="
			      "\"monospace\">Single</font></td>\n");
	    else
	      if (printopt0 == 'H')      printf("<TD>Single</TD>");
	      else if (printopt0 == 'g') fprintf(f,"E.Sngl");
	      else	                 fprintf(f,"Sngl");
	    nbExon = 0;
	  } else 
	    if (init) {
	      if (printopt0 == 'h')
		vhtml.push_back(" <td align=\"center\">\n  <font face="
				"\"monospace\">Initial</font></td>\n");
	      else
		if (printopt0 == 'H')      printf("<TD>Initial</TD>");
		else if (printopt0 == 'g') fprintf(f,"E.Init");
		else	                   fprintf(f,"Init");
	      if(!forward) nbExon = 0;
	    } else 
	      if (term) {
		if (printopt0 == 'h')
		  vhtml.push_back(" <td align=\"center\">\n  <font face="
				  "\"monospace\">Terminal</font></td>\n");
		else
		  if (printopt0 == 'H')       printf("<TD>Terminal</TD>");
		  else  if (printopt0 == 'g') fprintf(f,"E.Term");
		  else	                      fprintf(f,"Term");
		if(forward)  nbExon = 0;
	      }
	      else {
		if (printopt0 == 'h')
		  vhtml.push_back(" <td align=\"center\">\n  <font face="
				  "\"monospace\">Internal</font></td>\n");
		else		
		  if (printopt0 == 'H')       printf("<TD>Internal</TD>");
		  else  if (printopt0 == 'g') fprintf(f,"E.Intr");
		  else	                      fprintf(f,"Intr");
	      }
	  
	  if (printopt0 == 'g')
	    printf("\t%d\t%d\t.\t%c\t%d\n",
		   Lend,Rend,((forward) ? '+' : '-'),abs(PhaseAdapt(state))-1);
	  else
	    if (printopt0 == 'h')
	      vhtml.push_back(" <td align=\"right\">\n"
			      "  <font face=\"monospace\">" + to_string(Lend)
			      + "</font></td>\n"
			      " <td align=\"right\">\n"
			      "  <font face=\"monospace\">" + to_string(Rend)
			      + "</font></td>\n"
			      " <td align=\"center\">\n"
			      "  <font face=\"monospace\">."
			      "</font></td>\n"
			      " <td align=\"center\">\n"
			      "  <font face=\"monospace\">" + 
			      ((forward) ? '+' : '-') + "</font></td>\n"
			      " <td align=\"center\">\n"
			      "  <font face=\"monospace\">" +
			      to_string(abs(PhaseAdapt(state))-1) + 
			      "</font></td>\n</tr>\n");
	    else
	      if (printopt0 == 'H') {
		printf("</TD><TD>%c</TD><TD>%7d</TD><TD>%7d</TD>",
		       ((forward) ? '+' : '-'),Lend,Rend);
		printf("<TD>%4d</TD>", Rend-Lend+1);
		
		if (init)
		  fprintf(f,"<TD>%+2d</TD>", ((forward) ? 1: -1));
		else {
		  Phase = ((forward) ?
			   PhaseAdapt(stateBack-IntronF1) :
			   -PhaseAdapt(stateNext-IntronR1));
		  
		  if (abs(Phase) <= 3 && abs(Phase) >=1) fprintf(f,"<TD>%+2d</TD>",Phase);
		  else fprintf(f,"<TD>Unk.</TD>");
		}
		fprintf(f,"<TD>%+2d</TD>\n   </TR>\n",PhaseAdapt(state));
	      }
	      else {
		fprintf(f,"    %c    %7d %7d",((forward) ? '+':'-'),Lend,Rend);
		fprintf(f,"     %4d  ", Rend-Lend+1);
		
		if (init)
		  fprintf(f,"   %+2d", ((forward) ? 1: -1));
		else {
		  Phase = ((forward) ?
			   PhaseAdapt(stateBack-IntronF1) :
			   -PhaseAdapt(stateNext-IntronR1));
		  
		  if (abs(Phase) <= 3 && abs(Phase) >=1) fprintf(f,"   %+2d",Phase);
		  else fprintf(f," Unk.");
		}
		fprintf(f,"      %+2d",PhaseAdapt(state));
		fprintf(f," %7d %7d ", Don,Acc);
		fprintf(f,"  %3.0f.%-3.0f\n",100.0*(double)cons/(Rend-Lend+1),
			100.0*(double)incons/(Rend-Lend+1));
	      }
	}
	else 
	  if ((state >= UTR5F) && (state <= UTR3R)) {
	    if(printopt0 == 'h') {
	      lc++;
	      vhtml.push_back("<tr class=\"A" + to_string(lc%2) + "\">\n"+
			      " <td align=\"center\">\n"+
			      "  <font face=\"monospace\">" + position + "_" + 
			      to_string(nbGene) + "." + to_string(nbExon) +
			      "</font></td>\n");
	      vhtml.push_back(" <td align=\"center\">\n"
			      "  <font face=\"monospace\">");
	      vhtml.push_back(argv[0]);
	      vhtml.push_back("</font></td>\n");
	    }
	    else
	      if(printopt0 == 'g' || printopt0 == 'a')
		fprintf(f,nameformat,X->Name);
	      else
		if(printopt0 == 'H') {
		  lc++;
		  printf("   <TR class=A%d align=center>\n"
			 "    <TD>%d</TD><TD>%d</TD>",lc%2,nbGene,nbExon);
		}
		else {
		  fprintf(f,nameformat,X->Name);
		  fprintf(f,".%d.%d",nbGene,nbExon);
		}

	    if (printopt0 == 'g') printf("\tEuGene\t");
	    else if (printopt0 == 'a') printf(" ");
	    else if (printopt0 != 'H' && printopt0 != 'h') fprintf(f,"\t");
	    
	    switch (state) {
	    case UTR5F: // UTR5' F
	      if (printopt0 == 'h')
		vhtml.push_back(" <td align=\"center\">\n"
				"  <font face=\"monospace\">Utr5</font></td>\n"
				" <td align=\"right\">\n"
				"  <font face=\"monospace\">" +
				to_string(offset+posBack+1) + "</font></td>\n"
				" <td align=\"right\">\n"
				"  <font face=\"monospace\">" +
				to_string(offset+pos) + "</font></td>\n"
				" <td align=\"center\">\n"
				"  <font face=\"monospace\">.</font></td>\n"
				" <td align=\"center\">\n"
				"  <font face=\"monospace\">+</font></td>\n"
				" <td align=\"center\">\n"
				"  <font face=\"monospace\">.</font></td>\n"
				"</tr>\n");
	      else
	        if (printopt0 == 'g')
		  printf("UTR5\t%d\t%d\t.\t+\t.\n",offset+posBack+1,offset+pos);
		else if (printopt0 == 'H') printf("<TD>Utr5</TD><TD>+</TD>");
		else fprintf(f,"Utr5    +");
	      break;
	      
	    case UTR3F: // UTR 3' F
	      nbGene++;
	      if (printopt0 == 'h')
		vhtml.push_back(" <td align=\"center\">\n"
				"  <font face=\"monospace\">Utr3</font></td>\n"
				" <td align=\"right\">\n"
				"  <font face=\"monospace\">" +
				to_string(offset+posBack+1) + "</font></td>\n"
				" <td align=\"right\">\n"
				"  <font face=\"monospace\">" +
				to_string(offset+pos) + "</font></td>\n"
				" <td align=\"center\">\n"
				"  <font face=\"monospace\">.</font></td>\n"
				" <td align=\"center\">\n"
				"  <font face=\"monospace\">+</font></td>\n"
				" <td align=\"center\">\n"
				"  <font face=\"monospace\">.</font></td>\n"
				"</tr>\n");
	      else
		if (printopt0 == 'g')
		  printf("UTR3\t%d\t%d\t.\t+\t.\n",offset+posBack+1,offset+pos);
		else if (printopt0 == 'H') printf("<TD>Utr3</TD><TD>+</TD>");
		else fprintf(f,"Utr3    +");
	      break;
	      
	    case UTR5R: // UTR5' R
	      nbGene++;
	      if (printopt0 == 'h')
		vhtml.push_back(" <td align=\"center\">\n"
				"  <font face=\"monospace\">Utr5</font></td>\n"
				" <td align=\"right\">\n"
				"  <font face=\"monospace\">" +
				to_string(offset+posBack+1) + "</font></td>\n"
				" <td align=\"right\">\n"
				"  <font face=\"monospace\">" +
				to_string(offset+pos) + "</font></td>\n"
				" <td align=\"center\">\n"
				"  <font face=\"monospace\">.</font></td>\n"
				" <td align=\"center\">\n"
				"  <font face=\"monospace\">-</font></td>\n"
				" <td align=\"center\">\n"
				"  <font face=\"monospace\">.</font></td>\n"
				"</tr>\n");
	      else
		if (printopt0 == 'g')
		  printf("UTR5\t%d\t%d\t.\t-\t.\n",offset+posBack+1,offset+pos);
		else if (printopt0 == 'H') printf("<TD>Utr5</TD><TD>-</TD>");
		else fprintf(f,"Utr5    -");
	      break;
	      
	    case UTR3R:// UTR 3' R
	      if (printopt0 == 'h')
		vhtml.push_back(" <td align=\"center\">\n"
				"  <font face=\"monospace\">Utr3</font></td>\n"
				" <td align=\"right\">\n"
				"  <font face=\"monospace\">" +
				to_string(offset+posBack+1) + "</font></td>\n"
				" <td align=\"right\">\n"
				"  <font face=\"monospace\">" +
				to_string(offset+pos) + "</font></td>\n"
				" <td align=\"center\">\n"
				"  <font face=\"monospace\">.</font></td>\n"
				" <td align=\"center\">\n"
				"  <font face=\"monospace\">-</font></td>\n"
				" <td align=\"center\">\n"
				"  <font face=\"monospace\">.</font></td>\n"
				"</tr>\n");
	      else
		if (printopt0 == 'g')
		  printf("UTR3\t%d\t%d\t.\t-\t.\n",offset+posBack+1,offset+pos);
		else if (printopt0 == 'H') printf("<TD>Utr3</TD><TD>-</TD>");
		else fprintf(f,"Utr3    -");
	      break;
	    }
	    
	    if(printopt0 != 'g' && printopt0 != 'H' && printopt0 != 'h') {
	      fprintf(f,"    %7d %7d", offset+posBack+1, offset+pos);
	      fprintf(f,"     %4d  ",  pos - (posBack+1) +1);
	      fprintf(f,"   NA      NA");
	      fprintf(f,"      NA      NA ");
	      fprintf(f,"  %3.0f.%-3.0f\n",
		      100.0*(double)cons/((offset+pos)  -(offset+posBack+1)+1),
		      100.0*(double)incons/((offset+pos)-(offset+posBack+1)+1));
	 }
	    else if (printopt0 == 'H') {
	      printf("<TD>%7d</TD><TD>%7d</TD>", offset+posBack+1, offset+pos);
	      printf("<TD>%4d</TD>",  pos - (posBack+1) +1);
	      printf("<TD>NA</TD><TD>NA</TD>\n   </TR>\n");
	    }
	    if(-3 <= PhaseAdapt(stateNext) && PhaseAdapt(stateNext) <= -1)
	      nbExon = pred->nbExon(nbGene) + 1;
	  }
	if(pos == Data_Len)
	  break;
      }
   } else {
     int decalage;
     int line = 0;
     int state, posBack = 0, pos;
     
     decalage = ((argc == optind+1) ? offset : offset*(sequence-optind));
     
     for(i=pred->size()-1; i!=-1; i--) {
       if(i != pred->size()-1)
	 posBack   = pred->getPos(i+1);
       state = pred->getState(i);
       pos   = pred->getPos(i);
       
       if (state <= TermR3) {
	 fprintf(f,"%d %d ",decalage+posBack+1, decalage+pos);
	 line = 0;
       } else 
	 if(i != pred->size()-1 && state == InterGen) {
	   line = 1;
	   fprintf(f,"\n");
	 }
     }
     if (f==stdout) {
       if(line) fprintf(f,"\n");
       else fprintf(f,"\n\n");
     }
   }
  
  if (printopt0 == 'h'  &&  sequence+1 == argc) {   // -ph && last seq
    EndHTML();
  }
}

//-------------------------------------------------
// -ph print on stdout the begin of the HTML output
//-------------------------------------------------
void StartHTML() {
  printf("<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.01 Transitional//EN\">\n"
	 "<html>\n"
	 "  <head>\n"
	 "    <title>EuGène : Prediction</title>\n"
	 "    <link rel=\"SHORTCUT ICON\" href=\"WEB/Images/eg.jpg\">\n"
	 "    <link rel=\"STYLESHEET\" type=\"text/css\" "
	 "href=\"WEB/Style/eugene.css\">\n\n"
	 "    <script language=\"JavaScript1.2\" "
	 "src=\"WEB/Javascripts/coolmenus4.js\">\n"
	 "      /*****************************************************"
	 "*******************\n"
	 "      Copyright(c)2001 Thomas Brattli (webmaster@dhtmlcentral.com)\n"
	 "      DHTML coolMenus - Get it at coolmenus.dhtmlcentral.com\n"
	 "      Version 4.0_beta\n"
	 "      This script can be used freely as long as all copyright\n"
	 "      messages are intact.\n"
	 "      Extra info - Coolmenus reference/help - Extra links to help"
	 " files ****\n"
	 "      CSS help: "
	 "http://192.168.1.31/projects/coolmenus/reference.asp?m=37\n"
	 "      General: "
	 "http://coolmenus.dhtmlcentral.com/reference.asp?m=35\n"
	 "      Menu properties : "
	 "http://coolmenus.dhtmlcentral.com/properties.asp?m=47\n"
	 "      Level properties: "
	 "http://coolmenus.dhtmlcentral.com/properties.asp?m=48\n"
	 "      Background bar  : "
	 "http://coolmenus.dhtmlcentral.com/properties.asp?m=49\n"
	 "      Item properties : "
	 "http://coolmenus.dhtmlcentral.com/properties.asp?m=50\n"
	 "      *****************************************************"
	 "*******************/\n"
	 "    </script>\n"
	 "    <script language=\"JavaScript1.2\" "
	 "src=\"WEB/Javascripts/cm_addins.js\">\n"
	 "      /****************************************************"
	 "********************\n"
	 "      Copyright(c)2001 Thomas Brattli(webmaster@dhtmlcentral.com)\n"
	 "      Coolmenus add-in file for more advanced featuers..\n"
	 "      *****************************************************"
	 "*******************/\n"
	 "    </script>\n"
	 "    <script language=\"JavaScript1.2\" "
	 "src=\"WEB/Javascripts/diap.js\">\n"
	 "    </script>\n"
	 "  </head>\n\n"
	 "  <body leftmargin=\"0\" topmargin=\"0\" marginwidth=\"0\" "
	 "marginheight=\"0\"\n"
	 "	bgcolor=\"white\">\n"
	 "    <script language=\"JavaScript1.2\" "
	 "src=\"WEB/Javascripts/euG_menu.js\">\n"
	 "    </script>\n\n"
	 "    <table height=\"100%%\" width=\"100%%\" "
	 "CELLSPACING=\"0\" CELLPADDING=\"0\"\n"
	 "	   border=\"0\" rows=\"2\">\n"
	 "	<tr height=\"140\">\n"
	 "	  <td BACKGROUND=\"WEB/Images/top.jpg\" "
	 "colspan=\"2\" height=\"140\"\n"
	 "	      valign=\"top\">&nbsp;</td>\n"
	 "	</tr>\n"
	 "	<tr>\n"
         "	  <td width=\"145\" valign=\"top\" align=\"left\"\n"
	 "	      BACKGROUND=\"WEB/Images/left.jpg\">\n"
	 "	    <img src=\"WEB/Images/left.jpg\" "
	 "width=\"145\" height=\"10\"></td>\n"
	 "	  <td width=\"100%%\" valign=\"top\">\n\n"
	 "	    <!-- DEBUT PAGE... -->\n"
	 "	    <table CELLSPACING=\"0\" CELLPADDING=\"6\" "
	 "border=\"0\" width=\"100%%\">\n"
	 "		<tr>\n"
	 "                  <td colspan=\"2\"> </td>\n"
	 "                </tr>\n"
	 "                <tr>\n"
	 "                  <td colspan=\"2\"> "
	 "<img src=\"WEB/Images/euGpred_on.jpg\">\n"
	 "                    <br> <br> </td>\n"
	 "                </tr>\n"
	 "                <tr>\n"
	 "		   <td>\n"
	 "		    <form name=\"formname\">\n"
	 "		      <div align=\"center\">\n"
	 "			<center>\n"
	 "			  <table cellspacing=\"3\" cellpadding=\"5\" "
	 "border=\"0\" bgcolor=\"#c0dbe2\">\n"
	 "			      <tr>\n"
	 "				<td colspan=\"3\" align=\"center\" "
	 "bgcolor=\"white\">\n");
}

//-------------------------------------------------
// -ph print on stdout the begin of the HTML output
//-------------------------------------------------
void EndHTML() {
  for (int i = 0; i<(int)vhtml.size(); i++)
    std::cout << vhtml[i];
  printf("                                              </table>\n"
	 "					    </td>\n"
	 "					  </tr>\n"
	 "				      </table>\n"
	 "				    </div>\n"
	 "				  </div>\n"
	 "				</td>\n"
	 "				<td  bgcolor=\"#9bc7d0\"></td>\n"
	 "			      </tr>\n"
	 "			  </table>\n"
	 "			</center>\n"
	 "		      </div>\n"
	 "		    </form>\n"
	 "		  </td>\n"
	 "		</tr>\n"
	 "	    </table>\n"
	 "	  </td>\n"
	 "	</tr>\n"
	 "    </table>\n"
	 "    <!-- FIN PAGE -->\n"
	 "  </body>\n" 
	 "</html>\n");
}

//---------------------------
// Convertion int en string
//---------------------------
std::string to_string(int i)
{
  std::ostringstream os;  // création du flux
  os << i;           // insertion de l'int dans le flux
  return (os.str()); // extrait la valeur 
}

// -------------------------------------------------------------------------
// Verif coherence EST: calcul le nombre de nuc. coherents et
// incoherents avec les match est
// debut/fin/etat: debut et fin de la seq. dont l'etat est etat
// cons/incons: retour des valeurs
// -------------------------------------------------------------------------
void CheckConsistency(int debut, int fin, int etat, 
		      int * cons, int* incons, DNASeq *X, MasterSensor* ms)
{
  int i, con = 0, inc = 0;
  DATA dTMP;

  // les valeurs qui sont coherentes avec chaque etat
  const unsigned char Consistent[NbTracks] = {
    HitForward|MForward,    HitForward|MForward,    HitForward|MForward,
    HitReverse|MReverse,    HitReverse|MReverse,    HitReverse|MReverse,
    HitForward|MForward,    HitForward|MForward,    HitForward|MForward,
    HitReverse|MReverse,    HitReverse|MReverse,    HitReverse|MReverse,
    HitForward|MForward,    HitForward|MForward,    HitForward|MForward,
    HitReverse|MReverse,    HitReverse|MReverse,    HitReverse|MReverse,
    HitForward|MForward,    HitForward|MForward,    HitForward|MForward,
    HitReverse|MReverse,    HitReverse|MReverse,    HitReverse|MReverse,
    GapForward|MForward,    GapForward|MForward,    GapForward|MForward,
    GapReverse|MReverse,    GapReverse|MReverse,    GapReverse|MReverse,
    0,
    HitForward|MForward, HitForward|MForward|GapForward,
    HitReverse|MReverse, HitReverse|MReverse|GapReverse,
    GapForward|MForward, GapReverse|MReverse,
    GapForward|MForward, GapReverse|MReverse
  };
  
  const unsigned char MaskConsistent[NbTracks] = {
    Hit|Margin,    Hit|Margin,    Hit|Margin,    
    Hit|Margin,    Hit|Margin,    Hit|Margin,    
    Hit|Margin,    Hit|Margin,    Hit|Margin,    
    Hit|Margin,    Hit|Margin,    Hit|Margin,    
    Hit|Margin,    Hit|Margin,    Hit|Margin,    
    Hit|Margin,    Hit|Margin,    Hit|Margin,    
    Gap|Margin,    Gap|Margin,    Gap|Margin,
    Gap|Margin,    Gap|Margin,    Gap|Margin,
    0,
    Margin|Hit,    Margin|Hit|Gap,
    Margin|Hit,    Margin|Hit|Gap,
    Gap|Margin,    Gap|Margin,
    Gap|Margin,    Gap|Margin
  };

  if (debut == -1) debut = 0;
  
  for (i = debut; i <fin; i++) {
    
    ms->GetInfoSpAt(Type_Content, X, i, &dTMP);
    
    // y a t'il de l'info
    if (dTMP.EstMatch) {
      // y a t'il une info incoherente avec l'etat
      if (dTMP.EstMatch & ~MaskConsistent[etat]) 
	inc++;
      else if (dTMP.EstMatch & Consistent[etat]) 
	con++;
      else 
	inc++;
    }
  }
  *cons = con;
  *incons = inc;
}
