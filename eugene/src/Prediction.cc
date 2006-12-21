// ------------------------------------------------------------------
// Copyright (C) 2004 INRA <eugene@ossau.toulouse.inra.fr>
//
// This program is open source; you can redistribute it and/or modify
// it under the terms of the Artistic License (see LICENSE file).
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
//
// You should have received a copy of Artistic License along with
// this program; if not, please see http://www.opensource.org
//
// $Id$
// ------------------------------------------------------------------
// File:     Prediction.cc
// Contents: class Prediction
// ------------------------------------------------------------------

#include <iostream>
#include "Prediction.h"
#include "MSensor.h"
#include "DNASeq.h"

extern Parameters   PAR;

// -----------------------
//  Prediction: IsOriginal checks is the 1st gene of the current
//  alternative prediction does not appear either as a gene in the
//  optimal prediction or as an already predicted alternative gene
// -----------------------
bool Prediction :: IsOriginal(Prediction* optPred, std::vector <Prediction*>& altPreds, int seuil)
{
  Gene *thisGene = this->vGene[0];
  Gene *otherGene;
  int mismatch;

  // First check in the optimal prediction.
  int idx;

  for (idx = 0; idx < optPred->vGene.size(); idx++) {
    mismatch = thisGene->isDifferent(*(optPred->vGene[idx]), seuil);
    if ((mismatch >= 0) && (mismatch <= seuil))
      return false;
  }
  // Now alternative predictions. They contain just one gene
  
  for (idx = 0; idx < altPreds.size(); idx++) {
    mismatch = thisGene->isDifferent(*(altPreds[idx]->vGene[0]), seuil);
    if  ((mismatch >= 0) && (mismatch <= seuil))
      return false;
  }
  
  return true;
}

// -----------------------
//  ESTScan: fills ESTMatch with just Gap and Hit. Used to detect supported UTR
// -----------------------
void Prediction :: ESTScan()
{
    int i, NumEST;
    Hits *ThisEST = NULL, *AllEST = NULL;
    Block *ThisBlock = NULL;
    char tempname[FILENAME_MAX+1];
    FILE *fEST;
    
    ESTMatch = NULL;
    
    // If no EST data, do  nothing
    if (!(PAR.getI("Sensor.Est.use", 0, PAR.getI("EuGene.sloppy")))) return; 

    ESTMatch = new unsigned char[X->SeqLen+1];
    for (i = 0; i <= X->SeqLen; i++) ESTMatch[i] = 0;

    strcpy(tempname, PAR.getC("fstname"));
    strcat(tempname, ".est");
    fEST = FileOpen(NULL, tempname, "r", PAR.getI("EuGene.sloppy"));
    NumEST = 0;
    if (fEST)
      {
    	AllEST = AllEST->ReadFromFile(fEST, &NumEST, -1, 0,X->SeqLen);
        fclose(fEST);
	
	ThisEST = AllEST; // start from first one.
	while (ThisEST)
	  {
	    ThisBlock = ThisEST->Match;
	    while (ThisBlock) 
	      {
		//Mark Hit
		for (i = ThisBlock->Start; i <= ThisBlock->End; i++) ESTMatch[i] |= Hit;
		// A Gap ?
		if (ThisBlock->Prev != NULL) // Mark Gap
		  for (i = ThisBlock->Prev->End; i <= ThisBlock->Start; i++) ESTMatch[i] |= Gap;
		
		ThisBlock = ThisBlock->Next;
	      }
	    ThisEST = ThisEST->Next;
	  }
	delete AllEST;
      }
}

/*************************************************************
 **                       Feature Object                    **
 *************************************************************/
// ------------------------
//  Default constructor.
// ------------------------
Feature :: Feature ()
{
  number   = 0;
  state    = -1;
  start    = -1;
  end      = -1;
  strand   = '.';
  phase    = 9;
  framegff = 9;
  frame    = 9;
}

// ------------------------
//  Constructor.
// ------------------------
Feature :: Feature (signed char n, int s, int e, char st, int fr)
{
  number   = 0;
  state    = n;
  start    = s;
  end      = e;
  strand   = st;
  phase    = 9;
  framegff = 9;
  frame    = fr;
}

// ------------------------
//  Default destructor.
// ------------------------
Feature :: ~Feature ()
{
}


/*************************************************************
 **                         Gene Object                     **
 *************************************************************/
// ------------------------
//  Default constructor.
// ------------------------
Gene :: Gene ()
{
  complete = 0;
  isvariant = false;
  hasvariant = 0;
  clear();
}

// ------------------------
//  Default destructor.
// ------------------------
Gene :: ~Gene ()
{
  for (int i=0; i<nbFea(); i++) {
    delete  vFea[i];
    vFea[i] = NULL;
  }
  vFea.clear();
}
// ------------------------
//  Default cleaner.
// ------------------------
void Gene :: clear()
{
  exNumber   = exLength   = 0;
  inNumber   = inLength   = 0;
  utrLength  = mrnaLength = geneLength = 0;
  cdsStart   = cdsEnd     = trStart    = trEnd = -1;
}

// ------------------------
//  comparison: same CDS ?
// Assumes complete genes with similar stucture (number of exons, positions)
// returns the cumulated diffence in position or a negative number if different structure.
// ------------------------
int Gene :: isDifferent (const Gene& o, int threshold)
{
  int idxo = 0;
  int idxt = 0;
  int nDiff = 0;

  //first coding exon in each
  while (State2Status[o.vFea[idxo]->state] != 3) idxo++;
  while (State2Status[this->vFea[idxt]->state] != 3) idxt++;

  int shift = 0;

  while ((State2Status[o.vFea[idxo+shift]->state] != 2) && //not UTR
	 (State2Status[this->vFea[idxt]->state] != 2) && // not UTR
	 (idxo+shift < o.vFea.size()) &&
	 (idxt+shift < this->nbFea()))
    {
      if ((o.vFea[idxo+shift]->state == this->vFea[idxt+shift]->state) &&
	  (abs(o.vFea[idxo+shift]->start - this->vFea[idxt+shift]->start) <= threshold) &&
	  (abs(o.vFea[idxo+shift]->end - this->vFea[idxt+shift]->end)  <= threshold))
	{
	  nDiff += abs(o.vFea[idxo+shift]->start - this->vFea[idxt+shift]->start);
	  nDiff += abs(o.vFea[idxo+shift]->end - this->vFea[idxt+shift]->end);
	  shift++;
	}
      else return -1;
    }
  // we assume the gene structures satisfy the structural constraints
  // of complete genes. If the final exon is identical, it is a term
  // exon and the 2 CDS are identical.
  return nDiff;
}

// ------------------------
//  comparison: same CDS ?
// Assumes complete genes
// ------------------------
bool Gene :: operator== (const Gene& o)
{
  int idxo = 0;
  int idxt = 0;

  //first coding exon in each
  while (State2Status[o.vFea[idxo]->state] != 3) idxo++;
  while (State2Status[this->vFea[idxt]->state] != 3) idxt++;

  int shift = 0;

  while ((State2Status[o.vFea[idxo+shift]->state] != 2) && //not UTR
	 (State2Status[this->vFea[idxt]->state] != 2) && // not UTR
	 (idxo+shift < o.vFea.size()) &&
	 (idxt+shift < this->nbFea()))
    {
      if ((o.vFea[idxo+shift]->start != this->vFea[idxt+shift]->start) ||
	  (o.vFea[idxo+shift]->end != this->vFea[idxt+shift]->end) ||
	  (o.vFea[idxo+shift]->state != this->vFea[idxt+shift]->state))
	return false;
      shift++;
    }
  // we assume the gene structures staisfy the structural constraints
  // of complete genes. If the final exon is identical, it is a term
  // exon and the 2 CDS are identical.
  return true;
}

// ------------------------
//  Add Feature.
// ------------------------
void Gene :: AddFeature (signed char state, int start, int end)
{
  char strand = (PhaseAdapt(state) > 0) ? '+' : '-';
  int  frame  = PhaseAdapt(state);

  // If state == utr PhaseAdapt return 0 so
  if (state == UTR5F  ||  state == UTR3F) strand = '+';
  
  vFea.push_back(new Feature(state, start, end, strand, frame));

  // Complete : +1 pour un exon init|sngl et +2 pour term|sngl
  if (state <= TermR3) {
    if (state <= InitR3)      complete += 1;
    else if (state <= SnglR3) complete += 3;
    else if (state >= TermF1) complete += 2;
  }
}

// ------------------------
//  Update.
// ------------------------
void Gene :: Update ()
{
  int state, phase, framegff, lengff = 0;
  int stateBack = -1;
  int stateNext = -1;
  int forward   = (vFea[0]->strand == '+') ? 1 : 0;
  
  clear();
  
  // Compute Gene informations
  for (int i=0; i<nbFea(); i++) {
    state = vFea[i]->state;
    if (state <= TermR3) {
      exNumber++;
      exLength += vFea[i]->end - vFea[i]->start + 1;
      if(cdsStart==-1) cdsStart = vFea[i]->start - 1;
      cdsEnd = vFea[i]->end - 1;
    }
    else if (state <= IntronR2AG  ||
	     (state <= IntronU3R && state >= IntronU5F)) {
      inNumber++;
      inLength += vFea[i]->end - vFea[i]->start + 1;
    }
    else if (state <= UTR3R) {
      utrLength += vFea[i]->end - vFea[i]->start + 1;
	}
  }
  trStart = vFea[0]->start - 1;
  trEnd = vFea[nbFea()-1]->end - 1;
  if(vFea[nbFea()-1]->state < UTR5F) cdsEnd = trEnd;

  mrnaLength = exLength   + utrLength;
  geneLength = mrnaLength + inLength;
  if(trStart == -1) { cdsStart = trStart = 0; }
  if(trEnd   == -1) { cdsEnd   = trEnd   = 0; }
  
  // Compute number (feature)
  int nb = (forward) ? 1 : exNumber;
  for (int i=0; i<nbFea(); i++) {
    state = vFea[i]->state;
    if (state <= TermR3) {
      vFea[i]->number = nb;
      (forward) ? nb++ : nb--;
    }
  }

  // Compute phase (feature)
  for (int i=0; i<nbFea(); i++) {
    state = vFea[i]->state;
    if (i+1 < nbFea()) stateNext = vFea[i+1]->state;
    if (forward)
      ((stateBack == -1) ? phase=-1 : phase = PhaseAdapt(stateBack));
    else
      ((stateNext == -1) ? phase=-1 : phase = PhaseAdapt(stateNext));

    if (state <= SnglR3) vFea[i]->phase = (forward) ? 1: -1;
    else {
      if (abs(phase) <= 6 && abs(phase) >= 4) 
	vFea[i]->phase = (forward) ? phase-3 : phase+3;
    }
    stateBack = state;
  }

  // Compute framegff (feature)
  for (int i=0; i<nbFea(); i++) {
    // On ne peut calculer la framegff qd :
    //  - on a pas le debut d'un gene forward
    //  - on a pas le gene entier pour un reverse
    if ( forward  && complete==2  || !forward && complete!=3) break;

    state = vFea[i]->state;
    if (state > TermR3) continue;

    if (!forward && state >= TermR1) lengff = exLength;

    if (state <= SnglR3) {
      vFea[i]->framegff = 0;
      if (forward) lengff = vFea[i]->end - vFea[i]->start +1;
    }
    else {
      if (!forward) lengff -= vFea[i]->end - vFea[i]->start +1;
      vFea[i]->framegff = ( 3 - (lengff % 3) ) %3;
      if (forward)  lengff +=  vFea[i]->end - vFea[i]->start +1;
    }
  }
}
 
// ------------------------
//  PrintGene.
// ------------------------
void Gene :: PrintInfo (FILE *F, int nb, char *seqName)
{
  bool nofirst = false;
  int i, state;
  char comp[10];
  switch (complete) {
  case 1:
    strcpy(comp, "Cfrag");
    break;
  case 2:
    strcpy(comp, "Nfrag");
    break;
  case 3:
    strcpy(comp, "Full");
    break;
  default:
    strcpy(comp, "?");
  }

  // CDS
  fprintf(F,"%s.%d  \tEuGene_misc\tCDS\t%d\t%d\t%d\t%c\t.\t",
	  seqName, nb, cdsStart+1, cdsEnd+1, exLength, vFea[0]->strand);
  fprintf(F,"%s\t", comp);
  for (i=0; i<nbFea(); i++) {
    if(vFea[i]->state <= TermR3) {
      if (nofirst) fprintf(F,",");
      else         nofirst = true;
      fprintf(F,"%d..%d",vFea[i]->start,vFea[i]->end);
    }
  }
  fprintf(F,"\n");
  
  // Gene
  fprintf(F,"%s.%d  \tEuGene_misc\tGene\t%d\t%d\t%d\t%c\t.\t",
	  seqName, nb, trStart+1, trEnd+1,
	  geneLength, vFea[0]->strand);
  fprintf(F,"%s\t", comp);
  nofirst = false;
  for (i=0; i<nbFea(); i++) {
    state = vFea[i]->state;
    if(state <= TermR3 || (state >= UTR5F && state <= UTR3R)) {
      char sep = ',';
      // Si exon et etat d'avant utr OU
      // si utr  et etat d'avant exon
      if((state <= TermR3 && (i!=0 && vFea[i-1]->state >= UTR5F)) ||
	 (state >= UTR5F  && (i!=0 && vFea[i-1]->state <= TermR3)))
	sep = ':';
      if (nofirst) fprintf(F,"%c",sep);
      else         nofirst = true;
      fprintf(F,"%d..%d",vFea[i]->start,vFea[i]->end);
    }
  }
  fprintf(F,"\n");
}


/*************************************************************
 **                        Prediction                       **
 ************************************************************/
// ------------------------
//  Default constructor.
// ------------------------
Prediction :: Prediction ()
{
  clear();
}

// ------------------------
//  Default constructor.
// ------------------------
Prediction :: Prediction (std::vector <int> vPos,
			  std::vector <signed char> vState)
{
  clear();
  
  int i;
  int start = 1;
  int nbPos = vPos.size()-1;  // == vState.size()-1;

  for(i=0; i<(int)vPos.size(); i++) {
    if(i!=0) start = vPos[i-1] + 1;

    // Si state=intergenique OU start > seqlen on passe
    if(vState[i] == InterGen || start > vPos[vPos.size()-1]) continue;
    
    // Si premier item OU nvx gene avec ou sans utr
    // (=> etat=(UTR5F|UTR3R) && etat precedent est intergenique)
    if(i==0 ||
       ((vState[i] == UTR5F || vState[i] == UTR3R) && vState[i-1] == InterGen) ||
       (vState[i] <= TermR3 && vState[i-1] == InterGen)) {
      nbGene++;
      vGene.push_back( new Gene() );
    }
    vGene[nbGene-1]->AddFeature(vState[i], start, vPos[i]);
  }
}

// ------------------------
//  Default destructor.
// ------------------------
Prediction :: ~Prediction ()
{
  for(int i=0; i<nbGene; i++) {
    delete vGene[i];
    vGene[i] = NULL;
  }
  vGene.clear();
}
// ------------------------
//  Initializer
// ------------------------
void Prediction :: clear()
{
  MS = NULL;
  X  = NULL;
  ESTMatch = NULL;
  seqName[0]  = '\000';
  nbGene      = 0;
  optimalPath = 0;
}
// --------------------------
//  TrimAndUpdate: modifies the prediction to Trim UTR according to EST support
//  Also computes all statistics (length and intron numbers) 
//  and deletes too short genes
// --------------------------
void Prediction :: TrimAndUpdate(DNASeq* x)
{
  int  state,start,end;
	
  X = x;
	
  if (PAR.getI("Output.UTRtrim")) {
	 
    // Fill in ESTMatch data if any
    ESTScan();
    
    for(int i=0; i<nbGene; i++) 
      {
	// delete unsupported or trim
	std::vector <Feature*>::iterator featindex;
	std::vector <Feature*>::reverse_iterator rfeatindex; 	
  		
	for (featindex = vGene[i]->vFea.begin(); featindex != vGene[i]->vFea.end(); )
	  {
	    state = (*featindex)->state;
	    if (state <= InterGen) break; // Not UTR or UTRIntron
	    start = (*featindex)->start;
	    end = (*featindex)->end;
	    if (UTRCheckAndTrim(&start,&end,state)) 
	      {
		//  				printf("Triming %d-%d (state %d) to S%d\n",(*featindex)->start,(*featindex)->end,
		//  					(*featindex)->state,start);
		(*featindex)->start = start;
		break;
	      }
	    else 
	      {
		//  				printf("Erasing %d-%d (state %d)\n",(*featindex)->start,(*featindex)->end,(*featindex)->state);
		featindex = vGene[i]->vFea.erase(featindex);
	      }
	  }
  		
	// Same on right
	for (rfeatindex = vGene[i]->vFea.rbegin(); rfeatindex != vGene[i]->vFea.rend(); )
	  {
	    state = (*rfeatindex)->state;
	    if (state <= InterGen) break; // Not UTR or UTRIntron
	    start = (*rfeatindex)->start;
	    end = (*rfeatindex)->end;
	    if (UTRCheckAndTrim(&start,&end,state)) 
	      {
		//  				printf("Triming %d-%d (state %d) to E%d\n",(*featindex)->start,(*featindex)->end,
		//  					(*featindex)->state,end);
		(*rfeatindex)->end = end;
		break;
	      } 
	    else 
	      {
		//  				printf("Erasing %d-%d (state %d)\n",(*rfeatindex)->start,(*rfeatindex)->end,(*rfeatindex)->state);
		rfeatindex++;
		vGene[i]->vFea.pop_back();
	      }
	  }
      }
  }
  	
  if (ESTMatch) delete [] ESTMatch; // ESTScan() delete
  ESTMatch = NULL;
	
	
  // Gene update and deletion of short genes
  std::vector <Gene*>::iterator geneindex;
  int MinCDSLen = PAR.getI("Output.MinCDSLen");
  int gIdx = 0;

  for (geneindex = vGene.begin(); geneindex != vGene.end(); )
    {
      int empty_gene = ((*geneindex)->nbFea() < 1 ) ? 1: 0;
      if ( ! empty_gene )
	(*geneindex)->Update();
      if ((*geneindex)->exLength <= MinCDSLen ||  empty_gene ) {
	nbGene--;
	geneindex = vGene.erase(geneindex);
      } else 
	{
	  (*geneindex)->geneNumber = gIdx;
	  gIdx++;
	  geneindex++;
	}
    }
}

// --------------------------
//  Delete genes outside of range
// --------------------------
void Prediction :: DeleteOutOfRange(int s,int e)
{
  std::vector <Gene*>::iterator geneindex;

  for (geneindex = vGene.begin(); geneindex != vGene.end(); )
    {
      if (((*geneindex)->trEnd < s) | ((*geneindex)->trStart > e))
        {
          nbGene--;
          geneindex = vGene.erase(geneindex);
        } else geneindex++;
    }
}

// --------------------------
//  print prediction (master)
// --------------------------
void Prediction :: Print (DNASeq* x, MasterSensor *ms, FILE *OPTIM_OUT, const char append)
{
 
  MS = ms;
  FILE *OUT;
  char nameformat[20];
  char outputFormat[20];
  int  trunclen =      PAR.getI("Output.truncate");
  strcpy(outputFormat, PAR.getC("Output.format"));

  char *mode;
  std::ios_base::openmode ccmode;

  if (append) 
    {
      mode = "ab";
      ccmode = std::ofstream::out |std:: ofstream::app;
    }
  else 
    {
      mode = "wb";
      ccmode = std::ofstream::out | std::ofstream::trunc;
    }

  // SeqName
  if (trunclen) sprintf(nameformat,"%%%d.%ds",trunclen,trunclen);
  else strcpy(nameformat,"%s");
  sprintf(seqName, nameformat, X->Name);
  
  // If optimisation mode
  if (OPTIM_OUT != NULL) {
    PrintEgnL(OPTIM_OUT, seqName);
  }
  else{
    for (int i=0; i<strlen(outputFormat); i++) {
      char filename[FILENAME_MAX];
      strcpy(filename,PAR.getC("prefixName"));
      switch (outputFormat[i]) {
      case 'a':
	OUT = FileOpen(NULL, strcat(filename, ".egn.ara"), mode);
	PrintEgnL(OUT, seqName, 1);
	fclose(OUT);
	break;
      case 'd':
	OUT = FileOpen(NULL, strcat(filename, ".egn.debug"), mode);
	PrintEgnD(OUT);
	fclose(OUT);
	break;
      case 'g':
	OUT = FileOpen(NULL, strcat(filename,".gff1"), mode);
	PrintGff(OUT, seqName);
	fclose(OUT);
//SEB
  {
    std::string filename_gff3(PAR.getC("prefixName"));
    filename_gff3 += ".gff3";
    std::ofstream out(filename_gff3.c_str(),ccmode);
  //	, std::ios_base::binary);
    //le 1 est a verifier
    out << Gff3Line::header(X->Name, 1, X->SeqLen);
    PrintGff3(out, seqName);
    out.close();
  }
//SEB
	break;
      case 'h':
	OUT = FileOpen(NULL, strcat(filename, ".html"), mode);
	PrintHtml(OUT, seqName);
	fclose(OUT);
	break;
      case 'l':
	OUT = FileOpen(NULL, strcat(filename, ".egn"), mode);
	PrintEgnL(OUT, seqName);
	fclose(OUT);
	break;
      case 's':
	OUT = FileOpen(NULL, strcat(filename, ".egn.short"), mode);
	PrintEgnS(OUT);
	fclose(OUT);
	break;
      case 'o':
	fprintf(stderr,"\n    Seq         Type    S       Lend    Rend   Length  Phase   Frame      Ac      Do      Pr\n\n");
	PrintEgnL(stdout, seqName);
	fprintf(stdout, "\n");
	break;
      }
    }
  }
}
// ------------------------
//  print prediction (GFF)
// ------------------------
void Prediction :: PrintGff (FILE *OUT, char *seqName)
{
  int offset = PAR.getI("Output.offset");
  int stepid = PAR.getI("Output.stepid");
  int estopt = PAR.getI("Sensor.Est.use", 0, PAR.getI("EuGene.sloppy"));
  int state, start, end;
  int incons = 0, cons = 0;
  
  /* name source feature start end score strand frame */
  for(int i=0; i<nbGene; i++) {
  	
    for(int j=0; j<vGene[i]->nbFea(); j++) {
      state = vGene[i]->vFea[j]->state;
      // Print introns ?
      if( (PAR.getI("Output.intron") == 0) &&
	  (state >= IntronF1 && state <= InterGen) ||
	  (state >= IntronU5F) )
	continue;
      
      start = vGene[i]->vFea[j]->start;
      end   = vGene[i]->vFea[j]->end;
      if (estopt) CheckConsistency(start, end, state, &cons, &incons);
               
      fprintf(OUT, "%s.%d",seqName, ((vGene[i]->geneNumber)*stepid)+1);
      if (vGene[i]->isvariant) fprintf(OUT, "%c",'a'+(vGene[i]->hasvariant-1));
      fprintf(OUT, ".%d\tEuGene\t%s\t%d\t%d\t%.0f.%.0f\t%c\t", vGene[i]->vFea[j]->number,
	      State2GFFString(state), start+offset, end+offset,
	      100.0*(double)cons/(end-start+1), 100.0*(double)incons/(end-start+1),
	      vGene[i]->vFea[j]->strand);
      if(vGene[i]->vFea[j]->framegff == 9)
		fprintf(OUT, ".\n");
      else
		fprintf(OUT, "%d\n",abs(vGene[i]->vFea[j]->framegff));
      }
    }
  //ajoute une ligne vide en fin de fichier
   fprintf(OUT, "\n");
}

//SEB ------------------------
//  print prediction (GFF3)
// ------------------------
void Prediction :: PrintGff3 (std::ofstream& out, char *seqName)
{
  int offset = PAR.getI("Output.offset");
  int stepid = PAR.getI("Output.stepid");
  int estopt = PAR.getI("Sensor.Est.use", 0, PAR.getI("EuGene.sloppy"));
  int state, start, end;
  int incons = 0, cons = 0;
  
  /* name source feature start end score strand frame */
  for(int i=0; i<nbGene; i++)
  {
//  --  --  la ligne du gene  --  --
    Gff3Line gene_line(seqName);
    //la ligne decrivant le gene
    gene_line.setType(SOFA_GENE);
    gene_line.setStart(vGene[i]->trStart+1);
    gene_line.setEnd(vGene[i]->trEnd+1);

    std::string gene_name = 
      ((vGene[i]->isvariant) ? 
       to_string(seqName) + '.' + to_string((vGene[i]->geneNumber)*stepid +1)+(char)('a'+(vGene[i]->hasvariant-1)) :
       to_string(seqName) + '.' + to_string((vGene[i]->geneNumber)*stepid +1));

   //l'identifiant du gene = type_sofa:nom_eugene
    std::string gene_id = Sofa::getName(SOFA_GENE)
                        + ":" + gene_name;
    gene_line.setAttribute("ID=" + gene_id);
    gene_line.addAttribute("Name=" + gene_name);
    gene_line.addAttribute("length=" + to_string(vGene[i]->geneLength));
    //je pars du principe que la premiere structure me donne le brin du gene
    if (vGene[i]->nbFea() > 0)
      gene_line.setStrand(vGene[i]->vFea[0]->strand);
    gene_line.print(out);
//  --  --  fin de la ligne du gene  --  --

//  --  --  la ligne de l'ARNm  -- <=> basee sur le gene --
    //En fait je me contente d'ecraser ce que je ne veux plus
    gene_line.setType(SOFA_MRNA);
    std::string rna_id = Sofa::getName(SOFA_MRNA) + ":"+ gene_name;
    gene_line.setAttribute("ID=" + rna_id);
    gene_line.addAttribute("nb_exon=" + to_string(vGene[i]->exNumber));
//----- fin de la ligne du mRNA

    Gff3Line *pre_arn_line, *arn_line;
    std::vector<Gff3Line*> pre_arn_lines, arn_lines;
    std::string fea_name(gene_name + ".");
    std::string fea_id;
    //postulat : ici on a un seul ARNm par gene
    int taille_mRNA = 0;

    for(int j=0; j<vGene[i]->nbFea(); j++) {
      state = vGene[i]->vFea[j]->state;
      // Print introns ?
      if( (PAR.getI("Output.intron") == 0) &&
          (state >= IntronF1 && state <= InterGen) ||
          (state >= IntronU5F) )
        continue;

      start = vGene[i]->vFea[j]->start;
      end   = vGene[i]->vFea[j]->end;
      if (estopt) CheckConsistency(start, end, state, &cons, &incons);

      int type_sofa = getTypeSofa(state);
      //selon le type on renseigne une ligne existante ou nouvelle
      // d'arn, de cds ou des 2
      switch(type_sofa)
      {
        //cas d'une region intergenique, ou d'un intron
        //on ne le signale que dans le pre_arn
        case SOFA_INTERGEN  :
        case SOFA_INTRON    :
        case SO_UTR_INTRON  :
              pre_arn_line = fillGff3Line(type_sofa, start+offset, end+offset,
                                           vGene[i]->vFea[j]->strand,
                                           vGene[i]->vFea[j]->framegff);
              //les attributs
              setGff3Attributes(pre_arn_line, state, type_sofa,
                                 fea_name, j, gene_id);
            //on l'ajoute au vecteur
              pre_arn_lines.push_back(pre_arn_line);
              break;
        //cas d'un UTR 3' ou 5'
        //c'est [la fin d']un exon du pre_arn, et un utr de l'arn
        case SOFA_5_UTR :
        case SOFA_3_UTR :
              if((!pre_arn_lines.empty()) &&
                  (previousExonMustBeUpdated(pre_arn_lines.back(), start+offset)))
              //fusion des exons
                pre_arn_lines.back()->setEnd(end+offset);
              else
              {
                pre_arn_line = fillGff3Line(SOFA_EXON, start+offset, end+offset,
                                           vGene[i]->vFea[j]->strand,
                                           vGene[i]->vFea[j]->framegff);
              //attributs
                setGff3Attributes(pre_arn_line, state, SOFA_EXON,
                                  fea_name, j, gene_id);
                pre_arn_lines.push_back(pre_arn_line);
              }
              
              arn_line = fillGff3Line(type_sofa, start+offset, end+offset,
                                        vGene[i]->vFea[j]->strand,
                                        vGene[i]->vFea[j]->framegff);
              //attributs
              setGff3Attributes(arn_line, state, type_sofa,
                                fea_name, j, rna_id);
              arn_lines.push_back(arn_line);
              //je mets à jour la taille de l'ARNm
              taille_mRNA += (end-start+1);
              break;
        //cas d'un exon
        //il est peut etre la fin du precedent
        case SOFA_EXON :
//               fin = -99999;
              if((!pre_arn_lines.empty()) &&//c'est le premier
                (previousExonMustBeUpdated(pre_arn_lines.back(), start+offset)))
              //fusion des exons
                pre_arn_lines.back()->setEnd(end+offset);
              else
              {
                pre_arn_line = fillGff3Line(SOFA_EXON, start+offset, end+offset,
                                              vGene[i]->vFea[j]->strand,
                                              vGene[i]->vFea[j]->framegff);
                pre_arn_lines.push_back(pre_arn_line);
              }
              //les attributs
//               je le mets ici pour permettre la MAJ d'ontology term.
//               Pour bien faire il faudrait une fonction SetOntologyTerm()
//                 que l'on placerait au niveau de la MAJ
//                 et le setAttributes resterait dans le else
              //si je suis sur un single exon, le code de l'ontology doit etre SOFA_EXON
              if((SnglF1<=state) && (state <=SnglR3))
              //je passe nbTracks pour etre sur d'etre hors zone
                setGff3Attributes(pre_arn_lines.back(), NbTracks, type_sofa, fea_name,
                                  vGene[i]->vFea[j]->number, gene_id);
              else
                setGff3Attributes(pre_arn_lines.back(), state, type_sofa, fea_name,
                                  vGene[i]->vFea[j]->number, gene_id);

              arn_line = fillGff3Line(SOFA_CDS, start+offset, end+offset,
                                      vGene[i]->vFea[j]->strand,
                                      vGene[i]->vFea[j]->framegff);
              //les attributs
              setGff3Attributes(arn_line, state, SOFA_CDS,
                                fea_name, vGene[i]->vFea[j]->number, rna_id);
              arn_lines.push_back(arn_line);
              //je mets à jour la taille de l'ARNm
              taille_mRNA += (end-start+1);
              break;
        //par defaut je n'ajoute que dans le pre-arn
        default :
                pre_arn_line = fillGff3Line(type_sofa, start+offset, end+offset,
                                            vGene[i]->vFea[j]->strand,
                                            vGene[i]->vFea[j]->framegff);
                //les attributs
                setGff3Attributes(pre_arn_line, state, SOFA_CDS, fea_name,
                                  j, gene_id);
                pre_arn_lines.push_back(pre_arn_line);
      }//fin switch
      //si on vient d'ajouter une ligne dans les ARNs
      //Je lui ajoute le score
      if((type_sofa==SOFA_5_UTR) ||
          (type_sofa==SOFA_3_UTR) ||
          (type_sofa==SOFA_EXON))
      {
        double score = 100.0*(double)cons/(end-start+1);
        char buff[6];   snprintf(buff, 6, "%.1f", score);
        arn_lines.back()->addAttribute("est_cons=" + std::string(buff));
        score  = 100.0*(double)incons/(end-start+1);
        snprintf(buff, 6, "%.1f", score);
        arn_lines.back()->addAttribute("est_incons=" + to_string(buff));
      }
    }//fin des structures du gene
      //on affiche les lignes, j'en profite pour liberer les pointeurs
      std::vector<Gff3Line*>::iterator itr;
      //le pre-arn
      for(itr = pre_arn_lines.begin();
          itr != pre_arn_lines.end();
          ++itr)
      {
        (*itr)->print(out);
        delete (*itr);
      }
      //J'ajoute la taille du mRNA à sa ligne
      gene_line.addAttribute("length=" + to_string(taille_mRNA));
      //j'affiche la ligne du mRNA juste avant les UTRs et CDS
      gene_line.print(out);
      //l'arn
      for(itr = arn_lines.begin();
          itr != arn_lines.end();
          ++itr)
      {
        (*itr)->print(out);
        delete (*itr);
      }
      //fin du gene
      out << Gff3Line::endComplexFeature();
    }//fin des genes
}
//SEB
// ----------------------------
//  print prediction (egn long)
// ----------------------------
void Prediction :: PrintEgnL (FILE *OUT, char *seqName, int a)
{
  int offset = PAR.getI("Output.offset");
  int stepid = PAR.getI("Output.stepid");
  int estopt = PAR.getI("Sensor.Est.use", 0, PAR.getI("EuGene.sloppy"));
  int state, forward, start, end, don, acc;
  int incons = 0, cons = 0;

  /* Seq Type S Lend Rend Length Phase Frame Ac Do Pr */
  for(int i=0; i<nbGene; i++) {
  	
    forward = (vGene[i]->vFea[0]->strand == '+') ? 1 : 0;
    for(int j=0; j<vGene[i]->nbFea(); j++) {
      state = vGene[i]->vFea[j]->state;

      // Print introns ?
      if( (PAR.getI("Output.intron") == 0) &&
	  (state >= IntronF1 && state <= InterGen) ||
	  (state >= IntronU5F) )
	continue;

      start = vGene[i]->vFea[j]->start;
      end   = vGene[i]->vFea[j]->end;
      
     if (estopt) CheckConsistency(start, end, state, &cons, &incons);
 
      if (forward) {
		don = offset + start - 1;
		acc = offset + end   + 1;
      } else {
		acc = offset + start - 1;
		don = offset + end   + 1;
      }
    	  
      // araset ?
      if (a) fprintf(OUT, "%s ",seqName);
      else   
	{
	  fprintf(OUT, "%s.%d",seqName, ((vGene[i]->geneNumber)*stepid)+1);
	  if (vGene[i]->isvariant)
	    fprintf(OUT, "%c", 'a'+(vGene[i]->hasvariant-1));
	  fprintf(OUT, ".%d\t", vGene[i]->vFea[j]->number);
	}

      fprintf(OUT, "%s    %c    %7d %7d     %4d  ",
	     State2EGNString(state),
	     vGene[i]->vFea[j]->strand,
	     start+offset, end+offset, end - start +1);

      if(state >= UTR5F) {
		fprintf(OUT,"   NA      NA");
		fprintf(OUT,"      NA      NA ");
      }
      else 
      {
		if(vGene[i]->vFea[j]->phase == 9)
	  		fprintf(OUT," Unk.");
		else
	  		fprintf(OUT,"   %+2d",vGene[i]->vFea[j]->phase);
	
		if(vGene[i]->vFea[j]->frame == 0)
	  		fprintf(OUT,"      NA");
		else
	  		fprintf(OUT,"      %+2d",vGene[i]->vFea[j]->frame);
	
		fprintf(OUT," %7d %7d ", don, acc);
      }
      fprintf(OUT,"  %3.0f.%-3.0f\n",100.0*(double)cons/(end-start+1),
	      100.0*(double)incons/(end-start+1));
    }
  }
}
	

// -----------------------------
//  print prediction (egn short)
// -----------------------------
void Prediction :: PrintEgnS (FILE *OUT)
{
  int offset = PAR.getI("Output.offset");
  int line = 0;
  int state, start, end;
  
  for(int i=0; i<nbGene; i++) {
    if (i!=0)
      fprintf(OUT,"\n");
    for(int j=0; j<vGene[i]->nbFea(); j++) {
      state = vGene[i]->vFea[j]->state;
      if (state <= TermR3) {
	start = offset + vGene[i]->vFea[j]->start;
	end   = offset + vGene[i]->vFea[j]->end;
	fprintf(OUT,"%d %d ", start, end);
	line = 0;
      }
    }
  }
}

// -----------------------------
//  print prediction (egn debug)
// -----------------------------
void Prediction :: PrintEgnD (FILE *OUT)
{
  DATA Data;

  fprintf(OUT, "   pos nt  EF1   EF2   EF3   ER1   ER2   ER3    IF    IR    IG   U5F   U5R   U3F   U3R   IUF   IUR FWD tSta  tSto   Sta   Sto   Acc   Don   Ins   Del REV tSta  tSto   Sta   Sto   Acc   Don   Ins   Del noF tSta  tSto   Sta   Sto   Acc   Don   Ins   Del noR tSta  tSto   Sta   Sto   Acc   Don   Ins   Del\n");
  for(int i=0; i<X->SeqLen; i++) {
    MS->GetInfoAt   (X, i, &Data);      
    MS->PrintDataAt (X, i, &Data, OUT);
  }
  fprintf(OUT, "   pos nt  EF1   EF2   EF3   ER1   ER2   ER3    IF    IR    IG   U5F   U5R   U3F   U3R   IUF   IUR FWD tSta  tSto   Sta   Sto   Acc   Don   Ins   Del REV tSta  tSto   Sta   Sto   Acc   Don   Ins   Del noF tSta  tSto   Sta   Sto   Acc   Don   Ins   Del noR tSta  tSto   Sta   Sto   Acc   Don   Ins   Del\n");
}

// --------------------------
//  print prediction (Html)
// --------------------------
void Prediction :: PrintHtml (FILE *OUT, char *seqName)
{
  int state, start, end;
  int offset     = PAR.getI("Output.offset");
  int stepid     = PAR.getI("Output.stepid");
  char *html_dir = new char[FILENAME_MAX+1];
  strcpy(html_dir, PAR.getC("web_dir"));

  StartHTML(html_dir, OUT);
  OutputHTMLFileNames(1, html_dir, OUT);
  fprintf(OUT,
	  "				  </select>\n"
	  "				</td>\n"
	  "				<td width=\"100\" "
	  "align=\"center\"\n"
	  "				    bgcolor=\"#c0dbe2\">\n"
	  "				  <img src=\"%s"
	  "/Images/next.jpg\"\n"
	  "				       onclick=\"next();\"\n"
	  "				       title=\"Next\" "
	  "align=\"middle\">\n"
	  "				  &nbsp; &nbsp; &nbsp;\n"
	  "				  <img src=\"%s"
	  "/Images/last.jpg\"\n"
	  "				       onclick=\"last();\"\n"
	  "				       title=\"Jump to end\" "
	  "align=\"middle\">\n", html_dir, html_dir);
  fprintf(OUT,
	  "				</td>\n"
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
  fprintf(OUT,
	  "<tr align=\"center\" class=\"fonce\">\n"
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
  if (nbGene == 0)
    fprintf(OUT,
	    "<tr class=\"A0\">\n <td colspan=\"8\" "
	    "align=\"center\">\n  <font face=\"monospace\">\n  "
	    "<b>No exons/genes predicted in your submitted "
	    "sequence !</b></td>\n</tr>\n");

  int lc = 0;
  for(int i=0; i<nbGene; i++) {
    for(int j=0; j<vGene[i]->nbFea(); j++) {
      state = vGene[i]->vFea[j]->state;
      
      // Print introns ?
      if( (PAR.getI("Output.intron") == 0) &&
	  (state >= IntronF1 && state <= InterGen) ||
	  (state >= IntronU5F) )
	continue;
      
      start = offset + vGene[i]->vFea[j]->start;
      end   = offset + vGene[i]->vFea[j]->end;
      fprintf(OUT,
	      "<tr class=\"A%d\">\n"
	      " <td align=\"center\">\n"
	      "  <font face=\"monospace\">%s_%d.%d</font></td>\n"
	      "<td align=\"center\">\n"
	      "  <font face=\"monospace\">EuGene_%s</font></td>\n",
	      lc++%2, seqName, (i*stepid)+1, vGene[i]->vFea[j]->number,
	      PAR.getC("EuGene.organism"));
      fprintf(OUT,
	      " <td align=\"center\">\n  <font face="
	      "\"monospace\">%s</font></td>\n", State2GFFString(state));
      fprintf(OUT,
	      " <td align=\"right\">\n"
	      "  <font face=\"monospace\">%d</font></td>\n"
	      " <td align=\"right\">\n"
	      "  <font face=\"monospace\">%d</font></td>\n"
	      " <td align=\"center\">\n"
	      "  <font face=\"monospace\">."
	      "</font></td>\n"
	      " <td align=\"center\">\n"
	      "  <font face=\"monospace\">%c</font></td>\n"
	      " <td align=\"center\">\n"
	      "  <font face=\"monospace\">",
	      start, end, vGene[i]->vFea[j]->strand);

      if(vGene[i]->vFea[j]->framegff == 9)
	fprintf(OUT, ".");
      else
	fprintf(OUT, "%d",abs(vGene[i]->vFea[j]->framegff));
      fprintf(OUT, "</font></td>\n</tr>\n");
    }
  }
  
  EndHTML(OUT);
}

// -----------------------------------------
//  State2EGNString (convert state to char*)
// -----------------------------------------
char* Prediction :: State2EGNString (int state)
{
  if(state <= InitR3)                  return "Init";
  if(state <= SnglR3)                  return "Sngl";
  if(state <= IntrR3)                  return "Intr";
  if(state <= TermR3)                  return "Term";
  if(state <= IntronR2AG)              return "Intron";
  if(state == InterGen)                return "InterG";
  if(state == UTR5F || state == UTR5R) return "Utr5";
  if(state == UTR3F || state == UTR3R) return "Utr3";
  if(state >= IntronU5F)               return "Intron";
}

// -----------------------------------------
//  State2GFFString (convert state to char*)
// -----------------------------------------
char* Prediction :: State2GFFString (int state)
{
  if(state <= InitR3)                  return "E.Init";
  if(state <= SnglR3)                  return "E.Sngl";
  if(state <= IntrR3)                  return "E.Intr";
  if(state <= TermR3)                  return "E.Term";
  if(state <= IntronR2AG)              return "Intron";
  if(state == InterGen)                return "InterG";
  if(state == UTR5F || state == UTR5R) return "UTR5";
  if(state == UTR3F || state == UTR3R) return "UTR3";
  if(state >= IntronU5F)               return "Intron";
}
// ----------------------------------------------------------------
// If the feature checked is a UTR, then the UTR will be checked for 
// support (bool returned) and start/end positions may be trimmed to 
// match with the first experimental (EST) evidence. All EST evidence
// is used.
// ----------------------------------------------------------------
bool Prediction :: UTRCheckAndTrim(int* debut, int* fin, int etat)
{	
	if ((etat == UTR5F) || (etat == UTR3R))
		for (int k = *debut; k <= *fin; k++)
		  if (GetESTMatch(k) & Hit) { *debut = k+1; return true; }
				
	if ((etat == UTR3F) || (etat == UTR5R))
		for (int k = *fin; k >= *debut; k--)
		  if (GetESTMatch(k) & Hit) { *fin = k+1; return true; }
			
	return false;	
}
// ----------------------------------------------------------------
// Verif coherence EST: calcul le nombre de nuc. coherents et
// incoherents avec les match est debut/fin/etat: debut et fin de la
// seq. dont l'etat est etat cons/incons: retour des valeurs.  Seuls
// les EST informant ESTmatch (non filtrs) sont pris en compte.
// ----------------------------------------------------------------
void Prediction :: CheckConsistency(int debut, int fin, int etat, 
				    int* cons, int* incons)
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
    GapForward|MForward,    GapForward|MForward,    GapForward|MForward,
    GapReverse|MReverse,    GapReverse|MReverse,    GapReverse|MReverse,
    0,
    HitForward|MForward, HitForward|MForward,
    HitReverse|MReverse, HitReverse|MReverse,
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
    Hit|Margin,    Hit|Margin,    Hit|Margin,    
    Hit|Margin,    Hit|Margin,    Hit|Margin,    
    Gap|Margin,    Gap|Margin,    Gap|Margin,
    Gap|Margin,    Gap|Margin,    Gap|Margin,
    Gap|Margin,    Gap|Margin,    Gap|Margin,
    Gap|Margin,    Gap|Margin,    Gap|Margin,
    0,
    Hit|Margin,    Hit|Margin,
    Hit|Margin,    Hit|Margin,
    Gap|Margin,    Gap|Margin,
    Gap|Margin,    Gap|Margin
  };
  
  for (i=Max(debut-1,0); i<fin; i++) {
    
    MS->GetInfoSpAt(Type_Content, X, i, &dTMP);
    
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

//---------------------------------------------------
// PrintHtml : -ph print the begin of the HTML output
//---------------------------------------------------
void Prediction :: StartHTML (char* html_dir, FILE *OUT) {
  char *d = new char[MAX_LINE];  GetStrDate(d);
  
  fprintf(OUT,
	  "<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.01 Transitional//EN\">\n"
	  "<html>\n"
	  "  <head>\n"
	  "    <title>EuGne : Prediction</title>\n"
	  "    <link rel=\"STYLESHEET\" type=\"text/css\" "
	  "href=\"%s/Style/eugene.css\">\n\n"
	  "    <script language=\"JavaScript1.2\" "
	  "src=\"%s/Javascripts/diap.js\">\n"
	  "    </script>\n"
	  "  </head>\n\n"
	  "  <body leftmargin=\"0\" topmargin=\"0\" marginwidth=\"0\" "
	  "marginheight=\"0\"\n"
	  "	bgcolor=\"white\">\n"
	  "    <table height=\"100%%\" width=\"100%%\" "
	  "CELLSPACING=\"0\" CELLPADDING=\"0\"\n"
	  "	   border=\"0\" rows=\"2\">\n"
	  "	<tr height=\"140\">\n"
	  "	  <td BACKGROUND=\"%s/Images/top.jpg\" "
	  "colspan=\"2\" height=\"140\"\n"
	  "	      valign=\"top\">&nbsp;</td>\n"
	  "	</tr>\n"
	  "	<tr>\n"
	  "	  <td width=\"145\" valign=\"top\" align=\"left\"\n"
	  "	      BACKGROUND=\"%s/Images/left.jpg\">\n"
	  "	    <img src=\"%s/Images/left.jpg\" "
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
	  "<img src=\"%s/Images/euGpred_on.jpg\">\n",
	  html_dir, html_dir, html_dir,
	  html_dir, html_dir, html_dir);
  
  fprintf(OUT,
	  "<br><font size=-1> - "
	  "<a href=\"http://www.inra.fr/bia/T/EuGene/index.html\">EuGene</a> "
	  "%s for %s, %s -",
	  PAR.getC("EuGene.version"), PAR.getC("EuGene.organism"), d);

  if(rindex(PAR.getC("Output.format"), 'g') != NULL)
    fprintf(OUT," <a href=\"./%s.gff\">Gff output</a> -", PAR.getC("prefixName"));

  if(PAR.getI("Sensor.Est.use")    && PAR.getI("Est.PostProcess")  ||
     PAR.getI("Sensor.BlastX.use") && PAR.getI("BlastX.PostProcess"))
    fprintf(OUT," <a href=\"./%s.misc_info\">Evidences</a> -", PAR.getC("prefixName"));

  fprintf(OUT,
	  "<br></font></td>\n"
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
	  "bgcolor=\"white\">\n",
	  PAR.getC("EuGene.version"), PAR.getC("EuGene.organism"), d);
  
  delete [] d;
}

//-------------------------------------------------
// PrintHtml : -ph print the end of the HTML output
//-------------------------------------------------
void Prediction :: EndHTML(FILE *OUT) {
  fprintf(OUT,
	  "                                              </table>\n"
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

// ------------------------
//  Print gene info.
// ------------------------
void Prediction :: PrintGeneInfo (FILE* F)
{
  fprintf(F, "#=========================================#\n");
  fprintf(F, "#=           Gene informations           =#\n");
  fprintf(F, "#=========================================#\n");
  for(int i=0; i<nbGene; i++) {
    vGene[i]->PrintInfo(F, i+1, seqName);
  }
}

// ------------------------------------
//  Find the index of the first gene overlapping a segment.
//  -1 if not found. Assumes a sorted prediction.
// ------------------------------------
Gene *Prediction :: FindGene (int start, int end)
{
  for(int i=0; i<nbGene; i++) {
    if (vGene[i]->trStart > end) return NULL;
    if (vGene[i]->trEnd < start) continue;
    return vGene[i];
  }
  return NULL;
}

// ------------------------
//  plot a prediction.
// ------------------------
void Prediction :: PlotPred ()
{
  const int predWidth = 2;
  int state, start, end;
  int endBack = 0;

  for(int i=0; i<nbGene; i++) {
    for(int j=0; j<vGene[i]->nbFea(); j++) {
      state = vGene[i]->vFea[j]->state;
      start = vGene[i]->vFea[j]->start;
      end   = vGene[i]->vFea[j]->end;
      
      // Plot interg      
      for (int k=endBack+1; k<start; k++)
	PlotBarI(k, 0, 0.4, predWidth, 4);
      endBack = end;
      
      for (int k=start; k<end; k++)
	PlotBarI(k, State2Phase[state], 0.4, predWidth, 1);
    }
  }
}

// -------------------------------
//  get the state for a given pos.
// -------------------------------
char Prediction :: GetStateForPos(int pos)
{
  for(int i=0; i<nbGene; i++) {
    if( pos < vGene[i]->vFea[0]->start ) return InterGen;
    else
      if( pos <= vGene[i]->vFea[vGene[i]->nbFea()-1]->end )
	for(int j=0; j<vGene[i]->nbFea(); j++)
	  if(pos <= vGene[i]->vFea[j]->end)
	    return vGene[i]->vFea[j]->state;
  }
  return -1;
}

// ------------------------
//  isStart.
// ------------------------
char* Prediction :: IsStart(int p)
{
  int state  = GetStateForPos (p);
  int nState = GetStateForPos (p+1);
    
  if(1 <= State2Phase[nState] && State2Phase[nState] <= 3  &&  state >= InterGen)
    return "True";
  if(-3 <= State2Phase[state] && State2Phase[state] <= -1  && 
     (nState >= InterGen || p+1 >= vGene[0]->cdsEnd + 1))
    return "True";
  return "False";
}

// ------------------------
//  isStop.
// ------------------------
char* Prediction :: IsStop(int p)
{
  int state  = GetStateForPos (p);
  int nState = GetStateForPos (p+1);
  
  if(state != -1  &&  1 <= State2Phase[state] && State2Phase[state] <= 3  &&  
     (nState >= InterGen || p+1 >= vGene[0]->cdsEnd + 1))
    return "True";
  if(-3 <= State2Phase[nState] && State2Phase[nState] <= -1 &&  state >= InterGen)
    return "True";
  return "False";
}

// ------------------------
//  isDon.
// ------------------------
char* Prediction :: IsDon(int p)
{
  int state  = GetStateForPos (p);
  int nState = GetStateForPos (p+1);
  
  if(1 <= State2Phase[state] && State2Phase[state] <= 3 &&  nState == IntronF1)
    return "True";
  if(-3 <= State2Phase[nState] && State2Phase[nState] <= -1 &&  state == IntronR1)
    return "True";
  return "False";
}

// ------------------------
//  isAcc.
// ------------------------
char* Prediction :: IsAcc(int p)
{
  int pState = GetStateForPos (p);
  int state  = GetStateForPos (p+1);
  
  if(1 <= State2Phase[state] && State2Phase[state] <= 3 &&  pState == IntronF1)
    return "True";
  if(-3 <= State2Phase[pState] && State2Phase[pState] <= -1  &&  state == IntronR1)
    return "True";
  return "False";
}
 
// ------------------------
// IsState: Is a nucleotide in a given state ? 
// BE CAREFUL: 
// this method could be used only for a prediction with one complete gene.
//         the prediction could be the representation of an external gff annotation 
//         that could not specify the UTR. In this case, the state from 0 to the first exon
//         in the annotation (the first element of vPos and vState)
//         is InterGen and the state after the last exon 
//         in the annotation (the last element of vPos and vState)
//         is not set (getStateForPos returns -1).
//         To identify Start Reverse and Stop Forward, an other condition is used
//         pos == vPos[0] that compares the position with 
//                             in Forward the end of the last exon of the gene (vPos[0])
//                             in Reverse the start of the first exon of the gene (vPos[0])
// ------------------------
bool Prediction :: IsState (DATA::SigType sig_type, int pos, char strand)
{
  bool is_state = false;
  bool bad_strand = false;
  int state, nState, pState;

  if (sig_type == DATA::Start) {
    state  = GetStateForPos (pos);
    nState = GetStateForPos (pos+1);
    if (strand == '+') {
      if ((1 <= State2Phase[nState] && State2Phase[nState] <= 3) && state == InterGen) 
	is_state = true;
    } else if (strand == '-') {
      if ((-3 <= State2Phase[state] && State2Phase[state] <= -1)  &&  
	  (nState == UTR5R || pos == vGene[0]->cdsEnd + 1)) 
	is_state = true;
    } else bad_strand = true;

  } else if (sig_type == DATA::Stop) {
    state  = GetStateForPos (pos);
    nState = GetStateForPos (pos+1);
    if (strand == '+') {
      if( (1 <= State2Phase[state] && State2Phase[state] <= 3) &&  
	  (nState == UTR3F || pos == vGene[0]->cdsEnd + 1))
      is_state = true;
    } else if (strand == '-') {
      if( (-3 <= State2Phase[nState] && State2Phase[nState] <= -1)  &&  
	  (state == UTR3R || state == InterGen) ) is_state = true;
    } else bad_strand = true;

  } else if (sig_type == DATA::Acc) {
    pState = GetStateForPos (pos);
    state  = GetStateForPos (pos+1);
    if (strand == '+') {
      if ((1 <= State2Phase[state] && State2Phase[state]<= 3) &&  
	  (4 <= State2Frame[pState] && State2Frame[pState] <= 6)) is_state = true;
    } else if (strand == '-') {
      if( (-3 <= State2Phase[pState] && State2Phase[pState] <= -1)  &&  
	  (-6 <= State2Frame[state] && State2Frame[state] <= -4)) is_state = true;
    } else bad_strand = true;

  } else if (sig_type == DATA::Don) {
    state  = GetStateForPos (pos);
    nState = GetStateForPos (pos+1);
    if (strand == '+') {
      if ((1 <= State2Phase[state] && State2Phase[state] <= 3) &&  
	  (4 <= State2Frame[nState] && State2Frame[nState] <= 6)) is_state = true;
    } else if (strand == '-') {
      if( (-3 <= State2Phase[nState] && State2Phase[nState] <= -1) &&  
	  (-6 <= State2Frame[state] && State2Frame[state] <= -4)) is_state = true;
    } else bad_strand = true;

  } else 
    {std::cerr<<"ERROR: bad state "<<sig_type<<" given in argument in Prediction::IsState.\n"; exit(2);}

  if (bad_strand)
    {std::cerr<<"ERROR: bad strand "<<strand<<"  given in argument in Prediction::IsState.\n"; exit(2);}

  return is_state;
}


//SEB--  --  --  --  --  --  --
//  --  --  delocalisation pour eclaircir le code
Gff3Line*
  Prediction::fillGff3Line(int type_sofa, int start, int end,
                            char strand, int framegff)
{
  Gff3Line* line = new Gff3Line(seqName);
  line->setType(type_sofa);
  line->setStart(start);
  line->setEnd(end);
  line->setStrand(strand);
  if(framegff != 9)
    line->setPhase(abs(framegff));
  return line;
}

bool
  Prediction::previousExonMustBeUpdated(Gff3Line* line, int start)
{
  if(!line)
std::cerr << "bad previous !" << std::endl;
//     return false;
  if(line->getType()==SOFA_EXON)
    return((line->getEnd()+1) == start);
  return false;
}

void
  Prediction::setGff3Attributes(Gff3Line* line, int type_egn,
                                int type_sofa, std::string fea_name,
                                int j, std::string gene_id)
{
  std::string fea_id = Sofa::getName(type_sofa) + ":"
                      + fea_name + to_string(j);
  line->setAttribute("ID=" + fea_id);
  line->addAttribute("Parent="+gene_id);
  line->addAttribute("Ontology_term="
                    + Sofa::codeToString(getTypeSofa(type_egn, false)));
}
