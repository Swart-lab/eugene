/*****************************************************************************/
/*             Copyright (c) 2002 by INRA. All rights reserved.              */
/*                 Redistribution is not permitted without                   */
/*                 the express written permission of INRA.                   */
/*                     Mail : tschiex@toulouse.inra.fr                       */
/*---------------------------------------------------------------------------*/
/* File         : EuGeneTk/SensorPlugins/SMachine/Sensor.SMachine.cc         */
/* Description  : Sensor SpliceMachine                                       */
/* Authors      : P.Bardou, S.Foissac, M.J.Cros, A.Moisan, T.Schiex          */
/* History      : Mars 2004                                                  */
/*****************************************************************************/

#include "Sensor.SMachine.h"
#include"../../EuGene/System.h"

#define NORM(x,n) (((n)+(Max(-(n),x)))/(n))

/*************************************************************
 **                        SensorSpliceMachine              **
 *************************************************************/

extern Parameters PAR;

// ----------------------
//  Default constructor.
// ----------------------
SensorSMachine :: SensorSMachine (int n, DNASeq *X) : Sensor(n)
{
  char *seqname;
  char tempname[FILENAME_MAX+1];
  struct stat DirStat;
  int error = 0;

  type = Type_Acc|Type_Don|Type_Start;
  
  fprintf(stderr, "Probing SpliceMachine (splice sites)..........");  
  fflush(stderr);

  seqname = PAR.getC("fstname");
  strcpy(tempname,seqname);
  strcat(tempname,".SMachine");
  if (stat(tempname,&DirStat))
    error = mkdir(tempname,0755);
  
  if (error) {
    perror("SpliceMachine directory: ");
    exit(2);
  }

  strcpy(tempname,seqname);
  strcat(tempname,".SMachine/");
  strcat(tempname, BaseName(seqname));
  strcat(tempname,".splices");
  if (!ProbeFile(tempname)) SpliceMachine();
    
  ReadSMachineF(tempname, X->SeqLen);
  fprintf(stderr,"forward,");
  fflush(stderr);

  strcpy(tempname,seqname);
  strcat(tempname,".SMachine/");
  strcat(tempname, BaseName(seqname));
  strcat(tempname,".splicesR");
  if (!ProbeFile(tempname)) SpliceMachine();

  ReadSMachineR(tempname, X->SeqLen);
  fprintf(stderr," reverse done\n");
  
  CheckSplices(X,vPosAccF, vPosDonF, vPosAccR, vPosDonR);

  // vectors for reverse are put in the increasing order
  reverse(vPosAccR.begin(), vPosAccR.end()); 
  reverse(vValAccR.begin(), vValAccR.end());
  reverse(vPosDonR.begin(), vPosDonR.end());
  reverse(vValDonR.begin(), vValDonR.end());

  fprintf(stderr, "Probing SpliceMachine (starts)................");  
  fflush(stderr);

  strcpy(tempname,seqname);
  strcat(tempname,".SMachine/");
  strcat(tempname, BaseName(seqname));
  strcat(tempname,".starts");
  if (!ProbeFile(tempname)) SpliceMachine();

  ReadStartF(tempname, X->SeqLen);
  fprintf(stderr,"forward,");

  strcpy(tempname,seqname);
  strcat(tempname,".SMachine/");
  strcat(tempname, BaseName(seqname));
  strcat(tempname,".startsR");
  if (!ProbeFile(tempname)) SpliceMachine();

  ReadStartR(tempname, X->SeqLen);
  fprintf(stderr," reverse done\n");

  CheckStart(X,vPosF, vPosR);

  // vectors for reverse are put in the increasing order
  reverse(vPosR.begin(), vPosR.end()); 
  reverse(vValR.begin(), vValR.end()); 
}

// ----------------------
//  Default destructor.
// ----------------------
SensorSMachine :: ~SensorSMachine ()
{
  vPosAccF.clear();  vPosAccR.clear();
  vPosDonF.clear();  vPosDonR.clear();
  vValAccF.clear();  vValAccR.clear();
  vValDonF.clear();  vValDonR.clear();

  vPosF.clear();     vValF.clear();
  vPosR.clear();     vValR.clear();
}
// ----------------------
//  Init SMachine.
// ----------------------
void SensorSMachine :: Init (DNASeq *X)
{
  accB = PAR.getD("SMachine.accB*",GetNumber());
  accP = PAR.getD("SMachine.accP*",GetNumber());
  donB = PAR.getD("SMachine.donB*",GetNumber());
  donP = PAR.getD("SMachine.donP*",GetNumber());

  startP = PAR.getD("SMachine.startP*",GetNumber());
  startB = PAR.getD("SMachine.startB*",GetNumber());

  indexR = indexF = 0;
  iAccF = iDonF = iAccR = iDonR = 0;
  PositionGiveInfo = -1;

  if (PAR.getI("Output.graph")) Plot(X);
}
// -----------------------------
//  Read NetGene2 forward file.
// -----------------------------
void SensorSMachine :: ReadSMachineF(char name[FILENAME_MAX+1], int SeqLen)
{
  FILE *fp;
  char buf[FILENAME_MAX];
  char sacc[10],sdon[10];
  char altsacc[10],altsdon[10];
  int i, j,pos;
  
  if (!(fp = fopen(name, "r"))) {
    fprintf(stderr, "cannot open splice sites file %s\n", name);
    exit(2);
  }
  fgets(buf,FILENAME_MAX-1,fp);
  
  for (i = 0; i < SeqLen; i++) {
    // final value of netgene2 
    j = fscanf(fp,"%d %*s %s %s %*s %*s %*s %*s %s %s %*s %*s",
	       &pos,altsdon,altsacc,sdon,sacc);

    if ((j < 4) || (pos != i+1)) {
      fprintf(stderr, "Error in splice sites file %s, line %d\n", name, i+2);
      exit(2);
    }

    if (sdon[0] == '-') strcpy(sdon,altsdon);
    if (sacc[0] == '-') strcpy(sacc,altsacc);
    
    if( atof(sacc) != 0.0 ) {
      vPosAccF.push_back( i+1 );
      vValAccF.push_back( atof(sacc) );
    }
    if( atof(sdon) != 0.0 ) {
      vPosDonF.push_back( i );
      vValDonF.push_back( atof(sdon) );
    }
  }
  fclose(fp);
}

// -----------------------------
//  Read NetGene2 reverse file.
// -----------------------------
void SensorSMachine :: ReadSMachineR(char *name, int SeqLen)
{
  FILE *fp;
  char buf[FILENAME_MAX];
  char sacc[10],sdon[10];
  char altsacc[10],altsdon[10];
  int i, j, k,pos;
  
  if (!(fp = fopen(name, "r"))) {
    fprintf(stderr, "cannot open splice sites file %s\n", name);
    exit(2);
  }
  fgets(buf,FILENAME_MAX-1,fp);

  k = SeqLen;
  for (i = 0; i < SeqLen; i++) {
    // final value of netgene2 
    j = fscanf(fp,"%d %*s %s %s %*s %*s %*s %*s %s %s %*s %*s",
	       &pos,altsdon,altsacc,sdon,sacc);
    if ((j < 4) || (pos != i+1)) {
      fprintf(stderr, "Error in splice sites file %s, line %d\n", name, i+2);
      exit(2);
     }
    if (sdon[0] == '-') strcpy(sdon,altsdon);
    if (sacc[0] == '-') strcpy(sacc,altsacc);

    if( atof(sacc) != 0.0 ) {
      vPosAccR.push_back( k-1 );
      vValAccR.push_back( atof(sacc));
    }
    if( atof(sdon) != 0.0 ) {
      vPosDonR.push_back( k );
      vValDonR.push_back( atof(sdon));
    }
    k--;
  }
  fclose(fp);
}
// --------------------------
//  Read start forward file.
// --------------------------
void SensorSMachine :: ReadStartF (char *name, int Len)
{
  FILE *fp;
  int i,j = -1;
  double force;

  if (!(fp = fopen(name, "r"))) {
    fprintf(stderr, "cannot open SpliceMachine start file %s\n", name);
    exit(2);
  }
  
  while (1) {
    i = fscanf(fp,"%d %lf %*s\n", &j, &force);
    if (i == EOF) break;
    if (i < 2) {
      fprintf(stderr, "Error in SpliceMachine start file %s, position %d\n", name, j);
      exit(2);
    }
    vPosF.push_back( j-1 );
    vValF.push_back( force );
  }
  if (j == -1) fprintf(stderr,"WARNING: empty SpliceMachine start file !\n");
  fclose(fp);
}

// --------------------------
//  Read start reverse file.
// --------------------------
void SensorSMachine :: ReadStartR (char *name, int Len)
{
  FILE *fp;
  int i,j = -1;
  double force;

  if (!(fp = fopen(name, "r"))) {
    fprintf(stderr, "cannot open SpliceMachine start file %s\n", name);
    exit(2);
  }

  while (1) {
    i = fscanf(fp,"%d %lf %*s\n", &j, &force);
    if (i == EOF) break;
    if (i < 2) {
      fprintf(stderr, "Error in SpliceMachine start file %s, position %d\n", name, j);
      exit(2);
    }
    j = Len-j+2;
    vPosR.push_back( j-1 );
    vValR.push_back( force );
  }
  if (j == -1) fprintf(stderr,"WARNING: empty SpliceMachine start file !\n");
  fclose(fp);
}

// ------------------------
//  GiveInfo signal SMachine.   pow(atof(sdon),  donB)*donP
// ------------------------
void SensorSMachine :: GiveInfo (DNASeq *X, int pos, DATA *d)
{
  bool update = false;
  double f;

  if ( (PositionGiveInfo == -1) || (pos != PositionGiveInfo+1) ) update = true; // update indexes on vectors
  PositionGiveInfo = pos;

  // Accepteur Forward
  if(!vPosAccF.empty()) {
    if (update) 
      iAccF = lower_bound(vPosAccF.begin(), vPosAccF.end(), pos)-vPosAccF.begin();
    
    if((iAccF<(int)vPosAccF.size()) && (vPosAccF[iAccF] == pos)) {
      f = pow( vValAccF[iAccF], accB) * accP;
      d->sig[DATA::Acc].weight[Signal::Forward] += log(f);
      d->sig[DATA::Acc].weight[Signal::ForwardNo] += log(1.0-f);
      iAccF++;
    }
  }
  
  // Accepteur Reverse
  if (!vPosAccR.empty()) {
    if (update)
      iAccR = lower_bound(vPosAccR.begin(), vPosAccR.end(), pos)-vPosAccR.begin();

    if((iAccR<(int)vPosAccR.size()) && (vPosAccR[iAccR] == pos)) {
      f = pow( vValAccR[iAccR], accB) * accP;
      d->sig[DATA::Acc].weight[Signal::Reverse] += log(f);
      d->sig[DATA::Acc].weight[Signal::ReverseNo] += log(1.0-f);
      iAccR++;
    }
  }

  // Donneur Forward
  if(!vPosDonF.empty()) {
    if (update) 
      iDonF = lower_bound(vPosDonF.begin(), vPosDonF.end(), pos)-vPosDonF.begin();

    if ((iDonF<(int)vPosDonF.size()) && (vPosDonF[iDonF] == pos)) {
      f = pow( vValDonF[iDonF], donB) * donP;
      d->sig[DATA::Don].weight[Signal::Forward] += log(f);
      d->sig[DATA::Don].weight[Signal::ForwardNo] += log(1.0-f);
      iDonF++;
    }
  }
  
  // Donneur Reverse
  if(!vPosDonR.empty()) {
    if (update) 
      iDonR = lower_bound(vPosDonR.begin(), vPosDonR.end(), pos)-vPosDonR.begin();

    if((iDonR<(int)vPosDonR.size()) && (vPosDonR[iDonR] == pos)) {
      f = pow( vValDonR[iDonR], donB) * donP;
      d->sig[DATA::Don].weight[Signal::Reverse] += log(f);
      d->sig[DATA::Don].weight[Signal::ReverseNo] += log(1.0-f);
      iDonR++;
    }
  }

  // Start Forward
  if(!vPosF.empty()) {
    if (update) 
      indexF = lower_bound(vPosF.begin(), vPosF.end(), pos)-vPosF.begin();
    
    if((indexF<(int)vPosF.size()) && (vPosF[indexF] == pos)) {
      f = pow(vValF[indexF], startB)*(exp(-startP));
      d->sig[DATA::Start].weight[Signal::Forward] += log(f);
      d->sig[DATA::Start].weight[Signal::ForwardNo] += log(1.0-f);
      indexF++;
    }
  }
  
  // Start Reverse
  if (!vPosR.empty()) {
    if (update) 
      indexR = lower_bound(vPosR.begin(), vPosR.end(), pos)-vPosR.begin();
    
    if((indexR<(int)vPosR.size()) && (vPosR[indexR] == pos)) {
      f = pow(vValR[indexR], startB)*(exp(-startP));
      d->sig[DATA::Start].weight[Signal::Reverse] += log(f);
      d->sig[DATA::Start].weight[Signal::ReverseNo] += log(1.0-f);
      indexR++;
    }
  }
}

// ----------------------------
//  Plot Sensor information
// ----------------------------
void SensorSMachine :: Plot(DNASeq *X)
{
  double f;

  for (int i =0; i < (int)vPosAccF.size(); i++) {
    f = pow(vValAccF[i], accB) * accP;
    PlotBarF(vPosAccF[i], 4, 0.5, NORM(log(f),20.0), 4);
  }
  
  for (int i =0; i < (int)vPosDonF.size(); i++) {
    f = pow( vValDonF[i], donB) * donP;
    PlotBarF(vPosDonF[i], 4, 0.5, NORM(log(f),20.0), 11);
  }
  
  for (int i =0; i < (int)vPosAccR.size(); i++) {
    f = pow( vValAccR[i], accB) * accP;
    PlotBarF(vPosAccR[i], -4, 0.5, NORM(log(f),20.0), 4);
  }

  for (int i =0; i < (int)vPosDonR.size(); i++) {
    f = pow( vValDonR[i], donB) * donP;
    PlotBarF(vPosDonR[i], -4, 0.5, NORM(log(f),20.0), 11);
  }

  for (int i =0; i < (int)vPosF.size(); i++) {
    f = pow(vValF[i], startB)*(exp(-startP));
    PlotBarF(vPosF[i], (vPosF[i]%3)+1, 0.5, NORM(log(f),4.0), 2);
  }
  
  for (int i =0; i < (int)vPosR.size(); i++) {
    f = pow(vValR[i], startB)*(exp(-startP));
    PlotBarF(vPosR[i], -((X->SeqLen-vPosR[i])%3)-1, 0.5, NORM(log(f),4.0), 2);
  }

}

// ------------------
//  Post analyse
// ------------------
void SensorSMachine :: PostAnalyse(Prediction *pred)
{
}
// -----------------------
//  SpliceMachine launcher
// -----------------------
void SensorSMachine :: SpliceMachine()
{
  char *SMcmd=PAR.getC("SMachine.cmd");
  char s1[FILENAME_MAX+1];

  strcpy(s1, PAR.getC("fstname"));
  strcat(s1, ".SMachine/");
  strcat(s1, BaseName(PAR.getC("fstname")));

  printf("link %s %s\n",PAR.getC("fstname"),s1);
  link(PAR.getC("fstname"),s1);

  strcpy(s1,SMcmd);
  strcat(s1, " ");
  strcat(s1, PAR.getC("fstname"));
  strcat(s1, ".SMachine/");
  strcat(s1, BaseName(PAR.getC("fstname")));

  printf("%s\n",s1);
  system(s1);
  strcpy(s1, PAR.getC("fstname"));
  strcat(s1, ".SMachine/");
  strcat(s1, BaseName(PAR.getC("fstname")));

  printf("unlink %s\n",s1);
  unlink(s1);
  return;
}
