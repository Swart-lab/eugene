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
// File:     Sensor.AnnotaStruct.cc
// Contents: Sensor AnnotaStruct
// ------------------------------------------------------------------

#include "Sensor.AnnotaStruct.h"
#include <functional>
#include <stdio.h>
#include <stdlib.h>
bool BySigPos(const Signals *A, const Signals *B)
{ return (A->pos < B->pos); };

struct BySigPosF : public std::binary_function<Signals *,Signals *,bool> {
    bool operator()(Signals * x, Signals * y) { return BySigPos(x, y); }
};

bool ByConStart(const Contents *A, const Contents *B)
{ return (A->start < B->start); };

bool ByConEnd(const Contents *A, const Contents *B)
{ return (A->end < B->end); };

/*************************************************************
 **                     Signals object                      **
 *************************************************************/
// -------------------------
//  Default constructor.
// -------------------------
Signals :: Signals ()
{
  pos   = -1;
  type  = -1;
  edge  = -1;
  strcpy (score , "" );
}

// -------------------------
//  Default constructor.
// -------------------------
Signals :: Signals (int p, int t, int e, char s[20])
{
  pos   = p;
  type  = t;
  edge  = e;
  strcpy (score , s );
}

// -------------------------
//  Default destructor.
// -------------------------
Signals :: ~Signals () {}

// -------------------------
//  Print
//--------------------------
void Signals :: PrintS ()
{
  char t[7];
  char s = '+';

  switch (type) {
  case DATA::tStart : strcpy(t, "tStart"); break;
  case DATA::tStop  : strcpy(t, "tStop "); break;
  case DATA::Start  : strcpy(t, "Start "); break;
  case DATA::Stop   : strcpy(t, "Stop  "); break;
  case DATA::Don    : strcpy(t, "Acc   "); break;
  case DATA::Acc    : strcpy(t, "Don   "); break;
  case DATA::Ins    : strcpy(t, "Ins   "); break;
  case DATA::Del    : strcpy(t, "Del   "); break;
  }
  if (edge) s = '-';
  fprintf(stdout, "%d\t%s %c %s\n", pos, t, s, score);
}


/*************************************************************
 **                    Contents object                      **
 *************************************************************/
// -------------------------
//  Default constructor.
// -------------------------
Contents :: Contents ()
{
  start = -1;
  end   = -1;
  type  = -1;
  score = 0.0 ;
}

// -------------------------
//  Default constructor.
// -------------------------
Contents :: Contents (int sta, int e, int t, float s)
{
  start = sta;
  end   = e;
  type  = t;
  score = s;
}

// -------------------------
//  Default destructor.
// -------------------------
Contents :: ~Contents () {}

// -------------------------
//  Print
//--------------------------
void Contents :: PrintC ()
{
  char t[11];

  switch (type) {
  case DATA::ExonF1     : strcpy(t, "ExonF1    "); break;
  case DATA::ExonF2     : strcpy(t, "ExonF2    "); break;
  case DATA::ExonF3     : strcpy(t, "ExonF3    "); break;
  case DATA::ExonR1     : strcpy(t, "ExonR1    "); break;
  case DATA::ExonR2     : strcpy(t, "ExonR2    "); break;
  case DATA::ExonR3     : strcpy(t, "ExonR3    "); break;
  case DATA::IntronF    : strcpy(t, "IntronF   "); break;
  case DATA::IntronR    : strcpy(t, "IntronR   "); break;
  case DATA::InterG     : strcpy(t, "InterG    "); break;
  case DATA::UTR5F      : strcpy(t, "UTR5F     "); break;
  case DATA::UTR5R      : strcpy(t, "UTR5R     "); break;
  case DATA::UTR3F      : strcpy(t, "UTR3F     "); break;
  case DATA::UTR3R      : strcpy(t, "UTR3R     "); break;
  case DATA::IntronUTRF : strcpy(t, "IntronUTRF"); break;
  case DATA::IntronUTRR : strcpy(t, "IntronUTRR"); break;
  }
  fprintf(stdout, "%d\t%d\t%s %f\n", start, end, t, score);
}


/*************************************************************
 **                     SensorAnnotaStruct                  **
 *************************************************************/

extern Parameters PAR;

// ----------------------
// Default constructor.
// ----------------------
SensorAnnotaStruct :: SensorAnnotaStruct (int n, DNASeq *X) : Sensor(n)
{
  char tempname[FILENAME_MAX+1];
  char startRead[20] ,stopRead[20] ,accRead[20] ,donRead[20] ,tStartRead[20] ,tStopRead[20] ;
  char exonRead[20], intronRead[20], cdsRead[20];
  // all types are possible.
  type = Type_Any;

  fileExt   = PAR.getC("AnnotaStruct.FileExtension", GetNumber());
  
  strcpy(exonRead,PAR.getC("AnnotaStruct.Exon*",        GetNumber()));
  strcpy(intronRead,PAR.getC("AnnotaStruct.Intron*",    GetNumber()));
  strcpy(cdsRead,PAR.getC("AnnotaStruct.CDS*",         GetNumber()));

  if (exonRead[0] != 'i')
  {
  	exonPAR   = atof(exonRead);
	exonInline=0;
  }
  else
  {
	exonPAR   = 0;
	exonInline= 1;
  }
  if (intronRead[0] != 'i')
  {
  	intronPAR = atof(intronRead);
	intronInline = 0;
  }	
  else
  {
	intronPAR = 0;
  	intronInline = 1;
  }
  if (cdsRead[0] != 'i')
  {
	cdsPAR = atof(cdsRead);
     	cdsInline = 0;
  }
  else
  {
	cdsPAR   = 0;
	cdsInline = 1;
  }

  strcpy(startRead,PAR.getC("AnnotaStruct.Start*",     GetNumber()));
  strcpy(stopRead,PAR.getC("AnnotaStruct.Stop*",       GetNumber()));
  strcpy(accRead,PAR.getC("AnnotaStruct.Acc*",         GetNumber()));
  strcpy(donRead,PAR.getC("AnnotaStruct.Don*",         GetNumber()));
  strcpy(tStartRead,PAR.getC("AnnotaStruct.TrStart*",  GetNumber()));
  strcpy(tStopRead,PAR.getC("AnnotaStruct.TrStop*",    GetNumber()));

  strcpy(startPAR,  PAR.getC("AnnotaStruct.StartType",   GetNumber()));
  strcpy(stopPAR,   PAR.getC("AnnotaStruct.StopType",    GetNumber()));
  strcpy(accPAR,    PAR.getC("AnnotaStruct.AccType",     GetNumber()));
  strcpy(donPAR,    PAR.getC("AnnotaStruct.DonType",     GetNumber()));
  strcpy(tStartPAR, PAR.getC("AnnotaStruct.TrStartType", GetNumber()));
  strcpy(tStopPAR,  PAR.getC("AnnotaStruct.TrStopType",  GetNumber()));
  if (startRead[0] != 'i')
  {
 	strcat(startPAR,  startRead);
 	startInline=0;
  }
  else
  {
	startInline= 1;
  }
  if (stopRead[0] != 'i')
  {
  	strcat(stopPAR,  stopRead);
	stopInline = 0;
  }	
  else
  {
  	stopInline = 1;
  }
  if (accRead[0] != 'i')
  {
	strcat(accPAR, accRead);
     	accInline = 0;
  }
  else
  {
	accInline = 1;
  }
  if (donRead[0] != 'i')
  {
 	strcat(donPAR, donRead);
 	donInline=0;
  }
  else
  {
	donInline= 1;
  }
  if (tStartRead[0] != 'i')
  {
  	strcat(tStartPAR, tStartRead);
	tStartInline = 0;
  }	
  else
  {
  	tStartInline = 1;
  }
  if (tStopRead[0] != 'i')
  {
	strcat(tStopPAR,  tStopRead);
     	tStopInline = 0;
  }
  else
  {
	tStopInline = 1;
  }
  
  //fprintf(stderr, "Parameters : exonPAR %f  intronPAR %f cdsPAR %f startPAR %s stopPAR %s accPAR %s donPAR %s tStartPAR %s  tStopPAR %s\n",exonPAR, intronPAR ,cdsPAR, startPAR, stopPAR, accPAR, donPAR, tStartPAR,tStopPAR );
  //fflush(stderr);
  inputFormat_ = to_string(PAR.getC("AnnotaStruct.format", GetNumber(),1));
  
  fprintf(stderr, "Reading %s file....", fileExt);
  fflush(stderr);
  strcpy(tempname, PAR.getC("fstname"));
  strcat(tempname, ".");
  strcat(tempname, fileExt);
  if ( inputFormat_ == "GFF3" )
  {
    strcat(tempname,".gff3");
    char * filenameSoTerms = PAR.getC("Gff3.SoTerms",0,0);
    char * soTerms = new char[FILENAME_MAX+1];
    strcpy(soTerms , PAR.getC("eugene_dir"));
    strcat(soTerms , filenameSoTerms );

    GeneFeatureSet * geneFeatureSet = new GeneFeatureSet (tempname, soTerms);
    ReadAnnotaStructGff3(*geneFeatureSet, X->SeqLen);
    delete [] soTerms;
    delete geneFeatureSet;
  }
  else
  {
    ReadAnnotaStruct(tempname, X->SeqLen);
  }
  fprintf(stderr, "done\n");
  fflush(stderr);
  
  
  // Sort vSig by pos
  sort(vSig.begin(), vSig.end(), BySigPos);
  
  // Sort vCon by end and by start
  sort(vCon.begin(), vCon.end(), ByConEnd);
  stable_sort(vCon.begin(), vCon.end(), ByConStart);
  
  // Delete redundancy
  for(int i=0; i<(int)vSig.size(); i++) {
    if(i == 0) continue;
    if(vSig[i]->pos ==vSig[i-1]->pos  && vSig[i]->type ==vSig[i-1]->type &&
       vSig[i]->edge==vSig[i-1]->edge && vSig[i]->score==vSig[i-1]->score) {
      vSig.erase(i+vSig.begin());
      i--;
    }
  }

  for(int i=0; i<(int)vCon.size(); i++) {
    if(i == 0) continue;
    if(vCon[i]->start==vCon[i-1]->start && vCon[i]->end  ==vCon[i-1]->end &&
       vCon[i]->type ==vCon[i-1]->type  && vCon[i]->score==vCon[i-1]->score) {
      vCon.erase(i+vCon.begin());
      i--;
    }
  }

  // Patch because of non start/stop sites due to non full gene predicted
  // by fgeneshpasa.  If start stop is non canonical -> remove it
  for(int i=(int)vSig.size()-1; i>=0; i--) {
    if (vSig[i]->type == DATA::Start)
      if(vSig[i]->edge) {
	if (!X->IsEStart(vSig[i]->pos-1, -1)) { 
		vSig.erase(i+vSig.begin());
		fprintf(stderr,"Annotastruct: bad ATG on reverse strand\n");
	}
      }
      else {
	if (!X->IsEStart(vSig[i]->pos, 1))    { 
		vSig.erase(i+vSig.begin());
		fprintf(stderr,"Annotastruct: bad ATG on forward strand\n");
	}
      }
    if (vSig[i]->type ==  DATA::Stop)
      if(vSig[i]->edge) {
	if (!X->IsStop(vSig[i]->pos+2, -1)) { vSig.erase(i+vSig.begin()); }
      }
      else {
	if (!X->IsStop(vSig[i]->pos-3, 1))  { vSig.erase(i+vSig.begin()); }
      }
  }

  // FOR DEBUG
  std::vector <int> vPosAccF, vPosAccR;
  std::vector <int> vPosDonF, vPosDonR;
  std::vector <int> vPosStaF, vPosStaR;
  std::vector <int> vPosStoF, vPosStoR;
  for(int i=0; i<(int)vSig.size(); i++) {
    //vSig[i]->PrintS();
    switch (vSig[i]->type) {
    case DATA::Start : 
      if(vSig[i]->edge) vPosStaR.push_back(vSig[i]->pos);
      else vPosStaF.push_back(vSig[i]->pos);
      break;
    case DATA::Stop : 
      if(vSig[i]->edge) vPosStoR.push_back(vSig[i]->pos);
      else vPosStoF.push_back(vSig[i]->pos);
      break;
    case DATA::Acc :
      if(vSig[i]->edge) vPosAccR.push_back(vSig[i]->pos);
      else vPosAccF.push_back(vSig[i]->pos);
      break;
    case DATA::Don :
      if(vSig[i]->edge) vPosDonR.push_back(vSig[i]->pos);
    else vPosDonF.push_back(vSig[i]->pos);
      break;
    }
  }
  CheckStart  (X,vPosStaF, vPosStaR);
  CheckStop   (X,vPosStoF, vPosStoR);
  CheckSplices(X,vPosAccF, vPosDonF, vPosAccR, vPosDonR);
  vPosStaF.clear(); vPosStaR.clear();
  vPosStoF.clear(); vPosStoR.clear();
  vPosAccF.clear(); vPosDonF.clear(); vPosAccR.clear(); vPosDonR.clear();
  

}

// ----------------------
//  Default destructor.
// ----------------------
SensorAnnotaStruct :: ~SensorAnnotaStruct ()
{
  // Clear the data structures
  vSig.clear();
  vCon.clear();
}

// ----------------------
//  Init.
// ----------------------
void SensorAnnotaStruct :: Init (DNASeq *X)
{
  char exonRead[20], intronRead[20], cdsRead[20];
  char startRead[20] ,stopRead[20] ,accRead[20] ,donRead[20] ,tStartRead[20] ,tStopRead[20] ;
 
 strcpy(exonRead,PAR.getC("AnnotaStruct.Exon*",        GetNumber()));
  strcpy(intronRead,PAR.getC("AnnotaStruct.Intron*",    GetNumber()));
  strcpy(cdsRead,PAR.getC("AnnotaStruct.CDS*",         GetNumber()));

  if (exonRead[0] != 'i')
  {
  	exonPAR   = atof(exonRead);
	exonInline=0;
  }
  else
  {
	exonPAR   = 0;
	exonInline= 1;
  }
  if (intronRead[0] != 'i')
  {
  	intronPAR = atof(intronRead);
	intronInline = 0;
  }	
  else
  {
	intronPAR = 0;
  	intronInline = 1;
  }
  if (cdsRead[0] != 'i')
  {
	cdsPAR = atof(cdsRead);
     	cdsInline = 0;
  }
  else
  {
	cdsPAR   = 0;
	cdsInline = 1;
  }

  strcpy(startRead,PAR.getC("AnnotaStruct.Start*",     GetNumber()));
  strcpy(stopRead,PAR.getC("AnnotaStruct.Stop*",       GetNumber()));
  strcpy(accRead,PAR.getC("AnnotaStruct.Acc*",         GetNumber()));
  strcpy(donRead,PAR.getC("AnnotaStruct.Don*",         GetNumber()));
  strcpy(tStartRead,PAR.getC("AnnotaStruct.TrStart*",  GetNumber()));
  strcpy(tStopRead,PAR.getC("AnnotaStruct.TrStop*",    GetNumber()));

  strcpy(startPAR,  PAR.getC("AnnotaStruct.StartType",   GetNumber()));
  strcpy(stopPAR,   PAR.getC("AnnotaStruct.StopType",    GetNumber()));
  strcpy(accPAR,    PAR.getC("AnnotaStruct.AccType",     GetNumber()));
  strcpy(donPAR,    PAR.getC("AnnotaStruct.DonType",     GetNumber()));
  strcpy(tStartPAR, PAR.getC("AnnotaStruct.TrStartType", GetNumber()));
  strcpy(tStopPAR,  PAR.getC("AnnotaStruct.TrStopType",  GetNumber()));
  if (startRead[0] != 'i')
  {
 	strcat(startPAR,  startRead);
 	startInline=0;
  }
  else
  {
	startInline= 1;
  }
  if (stopRead[0] != 'i')
  {
  	strcat(stopPAR,  stopRead);
	stopInline = 0;
  }	
  else
  {
  	stopInline = 1;
  }
  if (accRead[0] != 'i')
  {
	strcat(accPAR, accRead);
     	accInline = 0;
  }
  else
  {
	accInline = 1;
  }
  if (donRead[0] != 'i')
  {
 	strcat(donPAR, donRead);
 	donInline=0;
  }
  else
  {
	donInline= 1;
  }
  if (tStartRead[0] != 'i')
  {
  	strcat(tStartPAR, tStartRead);
	tStartInline = 0;
  }	
  else
  {
  	tStartInline = 1;
  }
  if (tStopRead[0] != 'i')
  {
	strcat(tStopPAR,  tStopRead);
     	tStopInline = 0;
  }
  else
  {
	tStopInline = 1;
  }
  

  PosSigGiveInfo = -1;
  PosConGiveInfo = -1;
  iSig = 0;
  iCon = 0;
  
  if (PAR.getI("Output.graph")) Plot(X);
}

// --------------------------
//  Read start forward file.
// --------------------------
void SensorAnnotaStruct :: ReadAnnotaStruct(char name[FILENAME_MAX+1], int len)
{
  FILE  *fp;
  char  line[MAX_LINE];
  int   i;
  int   startC, endC, edge;   // C -> content
  int   startS, endS;         // S -> signal       need for reverse
  char  strand;
  char  phase[2];
  float scF;
  char  feature[20];
  char  scoreC[20];
  int   frame = -1;
 
  fp = FileOpen(NULL, name, "r", PAR.getI("EuGene.sloppy"));
  
  int j=0;
  while(fp  &&  fgets (line, MAX_LINE, fp) != NULL) {
    j++;
    if (line[0] != '#') {
      // GFF line : seqn source feature start end score strand phase
      i = sscanf(line, "%*s %*s %s %d %d %s %c %s",
		 feature, &startC, &endC, scoreC, &strand, phase);
      if (i < 6) {
	if (i==-1) {
	  if(j==1)
	    fprintf(stderr,"WARNING: empty AnnotaStruct file !...");
	}
	else {
	  fprintf(stderr, "Error in AnnotaStruct file %s, line %d.\n",name,j);
	  exit(2);
	}
      }
      else {
	
	/* Score ? */
	if (strcmp(scoreC, ".")) {
	  if (scoreC[0] == 'p' || scoreC[0] == 's') scF = atof(scoreC+1);
	  else                                      scF = atof(scoreC);
	}
	else scF = 0.0;

	/* Phase ? */
	if (strcmp(phase, ".")) {
	  if (strand == '+') { frame = (startC - 1)   % 3; }
	  else               { frame = (len   - endC) % 3; }
	  frame = (frame + atoi(phase)) % 3;
	}

	/* Strand ? */
	if      (strand == '+') {
	  edge = 0;
	  startS = startC;
	  endS   = endC;
	}
	else if (strand == '-') {
	  edge = 1;
	  startS = endC+1;
	  endS   = startC-1;
	}
	else {
	  fprintf(stderr, "WARNING: feature %s line %d strand unknown"
		  " => ignored.\n", feature, j);
	  continue;
	}
	startC--;
	endC--;
	
	/* Parse by feature */
	// Low level (contents OR signals)
	if     (!strcmp(feature, "trStart"))
	  vSig.push_back(new Signals(startS-1, DATA::tStart, edge, scoreC));
	else if(!strcmp(feature, "trStop"))
	  vSig.push_back(new Signals(startS,   DATA::tStop,  edge, scoreC));
	else if(!strcmp(feature, "start"))
	  vSig.push_back(new Signals(startS-1, DATA::Start,  edge, scoreC));
	else if(!strcmp(feature, "stop"))
	  vSig.push_back(new Signals(startS,   DATA::Stop,   edge, scoreC));
	else if(!strcmp(feature, "acc"))
	  vSig.push_back(new Signals(startS,   DATA::Acc,    edge, scoreC));
	else if(!strcmp(feature, "don"))
	  vSig.push_back(new Signals(startS-1, DATA::Don,    edge, scoreC));
	else if(!strcmp(feature, "ins"))
	  vSig.push_back(new Signals(startS,   DATA::Ins,    edge, scoreC));
	else if(!strcmp(feature, "del"))
	  vSig.push_back(new Signals(startS,   DATA::Del,    edge, scoreC));
	else if(!strcmp(feature, "exon"))
	  PushInCon(startC, endC, scF, strand, phase, frame);
	else if(!strcmp(feature, "intron"))
	  vCon.push_back(new Contents(startC, endC, DATA::IntronF+edge,scF));
	else if(!strcmp(feature, "utr5"))
	  vCon.push_back(new Contents(startC, endC, DATA::UTR5F+edge,  scF));
	else if(!strcmp(feature, "utr3"))
	  vCon.push_back(new Contents(startC, endC, DATA::UTR3F+edge,  scF));
	else if(!strcmp(feature, "utr")) {
	  vCon.push_back(new Contents(startC, endC, DATA::UTR5F+edge,  scF));
	  vCon.push_back(new Contents(startC, endC, DATA::UTR3F+edge,  scF));
	}
	else if(!strcmp(feature, "intronutr"))
	  vCon.push_back(new Contents(startC,endC,DATA::IntronUTRF+edge,scF));
	
	// High level (contents OR/AND signals)
	else if(!strcmp(feature, "E.Init")) 
	{
	  vSig.push_back(new Signals(startS-1, DATA::Start, edge, startPAR));
	  vSig.push_back(new Signals(endS,     DATA::Don,   edge, donPAR));
	  PushInCon(startC, endC, cdsPAR, strand, phase, frame);
	}
	else if(!strcmp(feature, "E.Intr")) {
	  vSig.push_back(new Signals(startS-1, DATA::Acc,   edge, accPAR));
	  vSig.push_back(new Signals(endS,     DATA::Don,   edge, donPAR));
	  PushInCon(startC, endC, cdsPAR, strand, phase, frame);
	}
	else if(!strcmp(feature, "E.Term")) {
	  vSig.push_back(new Signals(startS-1, DATA::Acc,   edge, accPAR));
	  vSig.push_back(new Signals(endS,     DATA::Stop,  edge, stopPAR));
       	  PushInCon(startC, endC, cdsPAR, strand, phase, frame);
	}
	else if(!strcmp(feature, "E.Sngl")) {
	  vSig.push_back(new Signals(startS-1, DATA::Start, edge, startPAR));
	  vSig.push_back(new Signals(endS,     DATA::Stop,  edge, stopPAR));
	  PushInCon(startC, endC, cdsPAR, strand, phase, frame);
	}
	else if(!strcmp(feature, "UTR5")) {
	  vSig.push_back(new Signals (startS-1,DATA::tStart,edge, tStartPAR));
	  vCon.push_back(new Contents(startC,endC,DATA::UTR5F+edge, cdsPAR));
	}
	else if(!strcmp(feature, "UTR3")) {
	  vSig.push_back(new Signals (endS,    DATA::tStop, edge, tStopPAR));
	  vCon.push_back(new Contents(startC,endC,DATA::UTR3F+edge, cdsPAR));
	}
	else if(!strcmp(feature, "UTR")) {
	  vSig.push_back(new Signals(startS-1, DATA::tStart, edge,tStartPAR));
	  vSig.push_back(new Signals(endS,     DATA::tStop,  edge,tStopPAR));
	  vCon.push_back(new Contents(startC,endC,DATA::UTR5F+edge, cdsPAR));
	  vCon.push_back(new Contents(startC,endC,DATA::UTR3F+edge, cdsPAR));
	}
	else if(!strcmp(feature, "Intron")) {
	  vCon.push_back(new Contents(startC,endC,DATA::IntronF+edge,intronPAR));
	}
	
	else if(!strcmp(feature, "E.Any")) {
	  PushInCon(startC, endC, exonPAR, strand, phase, frame);
	  vCon.push_back(new Contents(startC,endC,DATA::UTR5F+edge, exonPAR));
	  vCon.push_back(new Contents(startC,endC,DATA::UTR3F+edge, exonPAR));
	}
	else if(!strcmp(feature, "Intron.Any")) {
	  vSig.push_back(new Signals (startS-1, DATA::Don, edge, donPAR));
	  vSig.push_back(new Signals (endS,     DATA::Acc, edge, accPAR));
	  vCon.push_back(new Contents(startC,endC,DATA::IntronF+edge,intronPAR));
	  vCon.push_back(new Contents(startC,endC,DATA::IntronUTRF+edge,intronPAR));
	}
	else if(!strcmp(feature, "E.First")) {
	  PushInCon(startC, endC, exonPAR, strand, phase, frame);
	  vSig.push_back(new Signals (startS-1,DATA::tStart, edge, tStartPAR));
	  vCon.push_back(new Contents(startC,endC,DATA::UTR5F+edge, exonPAR));
	}
	else if(!strcmp(feature, "E.Last")) {
	  PushInCon(startC, endC, exonPAR, strand, phase, frame);
	  vSig.push_back(new Signals (endS,  DATA::tStop, edge, tStopPAR));
	  vCon.push_back(new Contents(startC,endC,DATA::UTR3F+edge, exonPAR));
	}
	else if(!strcmp(feature, "E.Extreme")) {
	  PushInCon(startC, endC, exonPAR, strand, phase, frame);
	  vSig.push_back(new Signals(startS-1, DATA::tStart, edge,tStartPAR));
	  vSig.push_back(new Signals(endS,     DATA::tStop,  edge,tStopPAR));
	  vCon.push_back(new Contents(startC,endC,DATA::UTR5F+edge, exonPAR));
	  vCon.push_back(new Contents(startC,endC,DATA::UTR3F+edge, exonPAR));
	}
	else
	  fprintf(stderr, "WARNING: feature %s line %d unknown => ignored.\n",
		  feature, j);
      }
    }
  }
  if (fp) fclose(fp);
}

//gff3
void SensorAnnotaStruct ::ReadAnnotaStructGff3(GeneFeatureSet & geneFeatureSet , int len)
{
  int   startC, endC, edge;   // C -> content
  int   startS, endS;         // S -> signal       need for reverse
  char  strand;
  char  phase[2];
  float scF;
  char  feature[50];
  string  onthology_term;
  string idSo;

  char  scoreC[20];
  int   frame = -1;
 
  int j=0;
  vector<GeneFeature *>::iterator it = geneFeatureSet.getIterator();
  int nbGeneFeature=geneFeatureSet.getNbFeature();
  int i=0;
  for ( i=0 ; i < nbGeneFeature ; i++, it++ )
  {
    j++;
    strcpy (scoreC, "s");
    strcpy (feature, (*it)->getType().c_str());
    startC = (*it)->getLocus()->getStart();
    endC = (*it)->getLocus()->getEnd();
    
    /* Score ? */
    scF = (*it)->getScore();
    //score by default for ins/del;
    strcat ( scoreC, to_string(scF).c_str() ); 
    
    strand = (*it)->getLocus()->getStrand();
    strcpy (phase, to_string((*it)->getPhase()).c_str() );
    onthology_term = (*it)->getAttributes()->getOntologyTerm();
    //recup code SO
    idSo=(*it)->getType();
    if ( idSo.find("SO:") == string::npos )
    {
      string tmp=GeneFeatureSet::soTerms_->getIdFromName(idSo);
      idSo=tmp;
    }
    /* Phase ? */
    if (strcmp(phase, "."))
    {
      if (strand == '+') { frame = (startC - 1) % 3; }
      else               { frame = (len - endC) % 3; }
      frame = (frame + atoi(phase)) % 3;
    }

    /* Strand ? */
    if (strand == '+') 
    {
      edge = 0;
      startS = startC;
      endS   = endC;
    }
    else if (strand == '-') 
    {
      edge = 1;
      startS = endC+1;
      endS   = startC-1;
    }
    else 
    {
      fprintf(stderr, "WARNING: feature %s line %d strand unknown"
	  " => ignored.\n", feature, j);
      continue;
    }
    startC--;
    endC--;

    //Recupere le score inline (scF) ou dans le fichier (cdsPAR)
    float scCds=0.0;
    (cdsInline == 1) ? scCds = scF : scCds = cdsPAR;
    
    float scIntron=0.0;
    (intronInline == 1) ? scIntron = scF : scIntron = intronPAR;
    
    float scExon=0.0;
    (exonInline == 1) ? scExon = scF : scExon = exonPAR;
    
    // Low level signals : use inline values
    if      ( idSo =="SO:0000315" ) //trStart : transcription_start_site 
    {
	GetScore (DATA::tStart,scF, scoreC);
	vSig.push_back(new Signals(startS-1, DATA::tStart, edge, scoreC));
    }
    else if ( idSo =="SO:0000616" ) //trStop : transcription_end_site
    {	
	GetScore (DATA::tStop,scF, scoreC);
       	vSig.push_back(new Signals(startS,   DATA::tStop,  edge, scoreC));
    }
    else if ( idSo =="SO:0000318" ) //start : start_codon
    {
	GetScore (DATA::Start,scF, scoreC);
        vSig.push_back(new Signals(startS-1, DATA::Start,  edge, scoreC));
    }
    else if ( idSo =="SO:0000319" ) //stop : stop_codon
    {
	GetScore (DATA::Stop,scF, scoreC);
        vSig.push_back(new Signals(startS,   DATA::Stop,   edge, scoreC));
    }
    else if ( idSo == "SO:0000164") //acc : three_prime_splice_site
    {
	GetScore (DATA::Acc,scF, scoreC);
        vSig.push_back(new Signals(startS,   DATA::Acc,    edge, scoreC));
    }
    else if ( idSo == "SO:0000163") //don : five_prime_splice_site
    {
	GetScore (DATA::Don,scF, scoreC);
        vSig.push_back(new Signals(startS-1, DATA::Don,    edge, scoreC));
    }
    else if ( idSo == "SO:0000366") //ins : insertion_site
      vSig.push_back(new Signals(startS,   DATA::Ins,    edge, scoreC));
    else if ( idSo == "SO:0000687") //del : deletion_junction
      vSig.push_back(new Signals(startS,   DATA::Del,    edge, scoreC));

    // High level (contents OR/AND signals)
    else if ( idSo == "SO:0000316" && onthology_term == "SO:0000196") 
      //CDS && five_prime_coding_exon_region == E.Init
    {
      GetScore (DATA::Start,scF, scoreC);
      vSig.push_back(new Signals(startS-1, DATA::Start, edge, scoreC));
      GetScore (DATA::Don,scF,scoreC );
      vSig.push_back(new Signals(endS,     DATA::Don,   edge, scoreC));
      PushInCon(startC, endC, scCds , strand, phase, frame);
    }
    else if ( idSo == "SO:0000316" && onthology_term == "SO:0000004") 
      //CDS && interior_coding_exon == E.Intr
    {
      GetScore (DATA::Acc,scF, scoreC);
      vSig.push_back(new Signals(startS-1, DATA::Acc,   edge, scoreC));
      GetScore (DATA::Don,scF, scoreC);
      vSig.push_back(new Signals(endS,     DATA::Don,   edge, scoreC));
      PushInCon(startC, endC, scCds, strand, phase, frame);
    }
    else if ( idSo == "SO:0000316" && onthology_term == "SO:0000197") //CDS && three_prime_coding_exon_region == E.Term
    {
      GetScore (DATA::Acc,scF, scoreC);
      vSig.push_back(new Signals(startS-1, DATA::Acc,   edge, scoreC));
      GetScore (DATA::Stop,scF, scoreC);
      vSig.push_back(new Signals(endS,     DATA::Stop,  edge, scoreC));
      PushInCon(startC, endC, scCds, strand, phase, frame);
    }
    else if ( idSo == "SO:0000316" && onthology_term == "SO:0005845") //CDS && single_exon == "E.Sngl"
    {
      GetScore (DATA::Start,scF, scoreC);
      vSig.push_back(new Signals(startS-1, DATA::Start, edge, scoreC));
      GetScore (DATA::Stop,scF, scoreC);
      vSig.push_back(new Signals(endS,     DATA::Stop,  edge, scoreC));
      PushInCon(startC, endC, scCds, strand, phase, frame);
    }
    else if ( idSo == "SO:0000316") //CDS
    {
      PushInCon(startC, endC, scCds, strand, phase, frame);
    }
    else if ( idSo == "SO:0000204" ) //five_prime_UTR
    {
      GetScore (DATA::tStart,scF,scoreC);
      vSig.push_back(new Signals (startS-1,DATA::tStart,edge, scoreC));

      vCon.push_back(new Contents(startC,endC,DATA::UTR5F+edge, scCds));
    }
    else if ( idSo == "SO:0000205" ) //three_prime_UTR
    {
      GetScore (DATA::tStop,scF, scoreC);
      vSig.push_back(new Signals (endS,    DATA::tStop, edge, scoreC));
      vCon.push_back(new Contents(startC,endC,DATA::UTR3F+edge, scCds));
    }
    else if ( idSo == "SO:0000203" ) //UTR
    {
      
      GetScore (DATA::tStart,scF, scoreC);
      vSig.push_back(new Signals(startS-1, DATA::tStart, edge, scoreC));
      GetScore (DATA::tStop,scF, scoreC);
      vSig.push_back(new Signals(endS,     DATA::tStop,  edge, scoreC));

      vCon.push_back(new Contents(startC,endC,DATA::UTR5F+edge, scCds));
      vCon.push_back(new Contents(startC,endC,DATA::UTR3F+edge, scCds));
    }
    else if ( idSo == "SO:0000188" && onthology_term == "SO:0000191") //intron not UTR !
    {
      vCon.push_back(new Contents(startC,endC,DATA::IntronF+edge,scIntron));
    }
    
    else if ( idSo == "SO:0000147" ) //E.Any 
    {
      PushInCon(startC, endC, scExon, strand, phase, frame);
      vCon.push_back(new Contents(startC,endC,DATA::UTR5F+edge, scExon));
      vCon.push_back(new Contents(startC,endC,DATA::UTR3F+edge, scExon));
    }
    else if ( idSo == "SO:0000188" ) // "Intron.Any"
    {
      GetScore (DATA::Don,scF, scoreC);
      vSig.push_back(new Signals (startS-1, DATA::Don, edge, scoreC));
      GetScore (DATA::Acc,scF, scoreC);
      vSig.push_back(new Signals (endS,     DATA::Acc, edge, scoreC));
      vCon.push_back(new Contents(startC,endC,DATA::IntronF+edge,scIntron));
      vCon.push_back(new Contents(startC,endC,DATA::IntronUTRF+edge,scIntron));
    }
    else if ( idSo == "SO:0000147" && onthology_term == "SO:0000200") //"E.First"
    {
      PushInCon(startC, endC, scExon, strand, phase, frame);
      vCon.push_back(new Contents(startC,endC,DATA::UTR5F+edge, scExon));
      GetScore (DATA::tStart,scF, scoreC);
      vSig.push_back(new Signals (startS-1,DATA::tStart, edge, scoreC));
    }
    else if ( idSo == "SO:0000147" && onthology_term == "SO:0000202") // "E.Last"
    {
      PushInCon(startC, endC, scExon, strand, phase, frame);
      vCon.push_back(new Contents(startC,endC,DATA::UTR3F+edge, scExon));
      GetScore (DATA::tStop,scF, scoreC);
      vSig.push_back(new Signals (endS,  DATA::tStop, edge, scoreC));
    }
    else if (idSo == "SO:0000147" && onthology_term.find("SO:0000202")!= string::npos && onthology_term.find("SO:0000200")!= string::npos) 
    {
      PushInCon(startC, endC, scExon, strand, phase, frame);
      vCon.push_back(new Contents(startC,endC,DATA::UTR5F+edge, scExon));
      vCon.push_back(new Contents(startC,endC,DATA::UTR3F+edge, scExon));
      
      GetScore (DATA::tStart,scF, scoreC);
      vSig.push_back(new Signals(startS-1, DATA::tStart, edge,scoreC));
      GetScore (DATA::tStop,scF, scoreC);
      vSig.push_back(new Signals(endS,     DATA::tStop,  edge,scoreC));
    }
    else
      fprintf(stderr, "WARNING: feature %s line %d unknown => ignored.\n",
	      feature, j);
    
    //fprintf(stderr, "END : feature %s line %d idSO : %s, Ontomlogy_term: %s.\n",
	   // feature,j, idSo.c_str(), onthology_term.c_str() );
    
  }
   
}

void SensorAnnotaStruct :: GetScore ( int type , float scF, char scoreC[20])
{
	if (type == DATA::tStart)
	{ 	//contient le type : p /s, si pas inline contiens le score du PAR
	    	strcpy (scoreC ,tStartPAR ); 
		if (tStartInline ==1)
		{
	 		strcat (scoreC , to_string(scF).c_str()); // ajout du score lu
		}
	}
	else if (type == DATA::tStop)
	{ 	
	    strcpy (scoreC ,tStopPAR ); 
		if (tStopInline ==1)
		{
			strcat (scoreC , to_string(scF).c_str()); // ajout du score lu
		}	

	}
	else if (type == DATA::Start)
	{ 	
		strcpy (scoreC ,startPAR );
		if (startInline ==1)
		{
			strcat (scoreC , to_string(scF).c_str()); 
		}	
	}
	else if (type == DATA::Stop)
	{ 	
		strcpy (scoreC ,stopPAR );
		if (stopInline ==1)
		{
			strcat (scoreC , to_string(scF).c_str()); 
		}	
	}
	else if (type == DATA::Acc)
	{ 	
		strcpy (scoreC ,accPAR );
		if (accInline ==1)
		{
			strcat (scoreC , to_string(scF).c_str()); 
		}	
	}
	else if (type == DATA::Don)
	{ 	
		strcpy (scoreC ,donPAR );
		if (donInline ==1)
		{
			strcat (scoreC , to_string(scF).c_str()); 
		}	
	}
	else if (type == DATA::Don)
	{ 	
		strcpy (scoreC ,donPAR );
		if (donInline ==1)
		{
			strcat (scoreC , to_string(scF).c_str()); 
		}	
	}

}
// ----------------
//  push_back con.
// ----------------
void SensorAnnotaStruct :: PushInCon(int d, int e, float sc,
				     char st, char p[2], int f)
{
  int k;
  if (st == '-') k = 3;
  else           k = 0;
  if (! strcmp(p, ".") ) {
    vCon.push_back(new Contents(d, e, DATA::ExonF1+k, sc));
    vCon.push_back(new Contents(d, e, DATA::ExonF2+k, sc));
    vCon.push_back(new Contents(d, e, DATA::ExonF3+k, sc));
  }
  else
    vCon.push_back(new Contents(d, e, f+k, sc));
}

// ------------------------
//  GiveInfo.
// ------------------------
void SensorAnnotaStruct :: GiveInfo (DNASeq *X, int pos, DATA *d)
{
  Signals TmpSig;
  bool  update = false;
  int   i, j, iConTMP;
  float k;
  
  /* Signals */
  // update indexes on vector ?
  if((PosSigGiveInfo == -1) || (pos != PosSigGiveInfo+1)) update = true;
  PosSigGiveInfo = pos;
  if (!vSig.empty()) {
    if(update) {
      TmpSig.pos = pos;
      iSig = lower_bound(vSig.begin(), vSig.end(),
			 &TmpSig, BySigPosF()) - vSig.begin();
    }
    while((iSig<(int)vSig.size()) && (vSig[iSig]->pos == pos)) {
      i = vSig[iSig]->type;
      j = vSig[iSig]->edge;
      if (vSig[iSig]->score[0] == 'p'  ||  vSig[iSig]->score[0] == 's')
	k = atof(vSig[iSig]->score + 1);
      else
	k = atof(vSig[iSig]->score);
      if (vSig[iSig]->score[0] == 'p') {     // Probability [0 1]
	d->sig[i].weight[j]   += log(k);
	d->sig[i].weight[j+2] += log(1-k);
      }
      else                                   // Score ]-oo +oo[
	d->sig[i].weight[j] += k;
      iSig++;
    }
  }

  /* Contents */
  // update indexes on vector ?
  if((PosConGiveInfo == -1) || (pos != PosConGiveInfo+1)) update = true;
  PosConGiveInfo = pos;
  
  if(!vCon.empty()) {
    if(update) {
      iCon = 0;
      while(iCon < (int)vCon.size()  &&  pos > vCon[iCon]->end) iCon++;
    }
    iConTMP = iCon;
    while((iConTMP < (int)vCon.size())   &&
	  (pos >= vCon[iConTMP]->start)  &&  (pos <= vCon[iConTMP]->end)) {
      d->contents[vCon[iConTMP]->type] += vCon[iConTMP]->score;
      iConTMP++;
    }
    while(iCon < (int)vCon.size()  &&  pos > vCon[iCon]->end) iCon++;
  }
}

// ----------------------------
//  Plot Sensor information.
// ----------------------------
void SensorAnnotaStruct :: Plot(DNASeq *X)
{
}

// ------------------
//  Post analyse.
// ------------------
void SensorAnnotaStruct :: PostAnalyse(Prediction *pred, FILE *MINFO)
{
}
