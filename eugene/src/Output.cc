#include "Output.h"

#include<iostream>

extern Parameters   PAR;

std::vector <std::string>  vhtml; // (-ph)

void Output (DNASeq *X, MasterSensor* ms, Prediction *pred, int sequence, int argc, char * argv[], FILE* f)
{
  int  i;
  DATA Data;
  int  Data_Len  = X->SeqLen;
  char printopt0 = PAR.getC("Output.format")[0]; 
  int  offset    = PAR.getI("Output.offset");
  int  estopt    = PAR.getI("Sensor.Est.use");

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
    printf("   pos nt  EF1   EF2   EF3   ER1   ER2   ER3    IF    IR    IG   U5F   U5R   U3F   U3R FW: tSta tSto  Sta  Sto  Acc  Don  Ins  Del REV: tSta tSto  Sta  Sto  Acc  Don  Ins  Del noFWD: tSta tSto  Sta  Sto  Acc  Don  Ins  Del noREV: tSta tSto  Sta  Sto  Acc  Don  Ins  Del\n");
    for(int i=0; i<Data_Len ; i++) {
      ms->GetInfoAt   (X, i, &Data);      ms->PrintDataAt (X, i, &Data);
    }
    printf("   pos nt  EF1   EF2   EF3   ER1   ER2   ER3    IF    IR    IG   U5F   U5R   U3F   U3R FW: tSta tSto  Sta  Sto  Acc  Don  Ins  Del REV: tSta tSto  Sta  Sto  Acc  Don  Ins  Del noFWD: tSta tSto  Sta  Sto  Acc  Don  Ins  Del noREV: tSta tSto  Sta  Sto  Acc  Don  Ins  Del\n");
    
  } else 
    if ((printopt0 == 'l') || (printopt0 == 'a') || (printopt0 == 'g') ||
	(printopt0 == 'h') || (printopt0 == 'H')) { 
      int nbGene  = 1;
      int nbExon  = 0;
      int cons = 0, incons = 0;
      int forward,init,term,Lend,Rend,Phase;
      int Don,Acc;
      char seqn[6]   = "";
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
      
      // position = strstr(argv[sequence],"/seq");
      if (f==stdout) {
	position = BaseName(argv[sequence]);
	if (position  == NULL)
	  strcpy(seqn,"     ");
	else {
	  // on enleve l'extension (.fasta)
	  if (char * suffix = rindex(position,'.')) *suffix = 0;
	  strncpy(seqn,position,5);
	  if(strlen(seqn) < 5 && printopt0 != 'g')
	    for(i=strlen(seqn); i<5; i++)
	      strcat(seqn, " ");
	  seqn[5] = '\0';
	}
      } else 
	  strcpy(seqn,"OPTIM");
      
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
	if(i != 0)
	  stateNext = pred->getState(i-1);
	
	if(pos == 0)
	  continue;
	
	if (estopt)
	  CheckConsistency(posBack, pos, state, &cons, &incons, X, ms);
	
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
	      printf("%s",seqn);
	    else
	      if(printopt0 == 'H') {
		lc++;
		printf("   <TR class=A%d align=center>\n"
		       "    <TD>%d</TD><TD>%d</TD>",lc%2,nbGene,nbExon);
	      }
	      else
		fprintf(f,"%s.%d.%d.%d",seqn,sequence-optind+1,nbGene,nbExon);
	  
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
	      if (printopt0 == 'H') printf("<TD>Single</TD>");
	      else fprintf(f,"Sngl");
	    nbExon = 0;
	  } else 
	    if (init) {
	      if (printopt0 == 'h')
		vhtml.push_back(" <td align=\"center\">\n  <font face="
				"\"monospace\">Initial</font></td>\n");
	      else
		if (printopt0 == 'H') printf("<TD>Initial</TD>");
		else fprintf(f,"Init");
	      if(!forward) nbExon = 0;
	    } else 
	      if (term) {
		if (printopt0 == 'h')
		  vhtml.push_back(" <td align=\"center\">\n  <font face="
				  "\"monospace\">Terminal</font></td>\n");
		else
		  if (printopt0 == 'H') printf("<TD>Terminal</TD>");
		  else fprintf(f,"Term");
		if(forward)  nbExon = 0;
	      }
	      else {
		if (printopt0 == 'h')
		  vhtml.push_back(" <td align=\"center\">\n  <font face="
				  "\"monospace\">Internal</font></td>\n");
		else		
		  if (printopt0 == 'H') printf("<TD>Internal</TD>");
		  else fprintf (f,"Intr");
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
			   PhaseAdapt(stateBack-6) :
			   -PhaseAdapt(stateNext-9));
		  
		  if (abs(Phase) <= 3) fprintf(f,"<TD>%+2d</TD>",Phase);
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
			   PhaseAdapt(stateBack-6) :
			   -PhaseAdapt(stateNext-9));
		  
		  if (abs(Phase) <= 3) fprintf(f,"   %+2d",Phase);
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
		printf("%s",seqn);
	      else
		if(printopt0 == 'H') {
		  lc++;
		  printf("   <TR class=A%d align=center>\n"
			 "    <TD>%d</TD><TD>%d</TD>",lc%2,nbGene,nbExon);
		}
		else
		  fprintf(f,"%s.%d.%d.%d",seqn,sequence-optind+1,nbGene,nbExon);
	    
	    if (printopt0 == 'g') printf("\tEuGene\t");
	    else if (printopt0 == 'a') printf(" ");
	    else if (printopt0 != 'H' && printopt0 != 'h') fprintf(f,"\t");
	    
	    switch (state) {
	    case 13: // UTR5' F
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
		  printf("Utr5\t%d\t%d\t.\t+\t.\n",offset+posBack+1,offset+pos);
		else if (printopt0 == 'H') printf("<TD>Utr5</TD><TD>+</TD>");
		else fprintf(f,"Utr5    +");
	      break;
	      
	    case 14: // UTR 3' F
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
		  printf("Utr3\t%d\t%d\t.\t+\t.\n",offset+posBack+1,offset+pos);
		else if (printopt0 == 'H') printf("<TD>Utr3</TD><TD>+</TD>");
		else fprintf(f,"Utr3    +");
	      break;
	      
	    case 15: // UTR5' R
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
		  printf("Utr5\t%d\t%d\t.\t-\t.\n",offset+posBack+1,offset+pos);
		else if (printopt0 == 'H') printf("<TD>Utr5</TD><TD>-</TD>");
		else fprintf(f,"Utr5    -");
	      break;
	      
	    case 16:// UTR 3' R
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
		  printf("Utr3\t%d\t%d\t.\t-\t.\n",offset+posBack+1,offset+pos);
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
	    if(stateNext >= ExonR1 && stateNext <= ExonR3)
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
       
       if (state <= ExonR3) {
	 fprintf(f,"%d %d ",decalage+posBack+1, decalage+pos);
	 line = 0;
       } else 
	 if(i != pred->size()-1 && state == InterGen5 || state == InterGen3) {
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
		      int * cons, int* incons, DNASeq *X, MasterSensor* ms)
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
    
    ms->GetInfoSpAt(Type_Content, X, i, &dTMP);
    
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
