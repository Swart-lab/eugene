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

#define NORM(x,n) Max(3.0,(x)+(n))/(n)

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

  type = Type_Acc|Type_Don|Type_Start;
  
  fprintf(stderr, "Probing SpliceMachine (splice sites)..........");  
  fflush(stderr);

  seqname = PAR.getC("fstname");
  strcpy(tempname,seqname);
  strcat(tempname,".spliceMAD");
  if (!ProbeFile(tempname)) SpliceMachine();
    
  ReadSMachineSplices(tempname, X->SeqLen);
  fprintf(stderr,"  done\n");
  
  CheckSplices(X,vPosAccF, vPosDonF, vPosAccR, vPosDonR);

  // vectors for reverse are put in the increasing order
  reverse(vPosAccR.begin(), vPosAccR.end()); 
  reverse(vValAccR.begin(), vValAccR.end());
  reverse(vPosDonR.begin(), vPosDonR.end());
  reverse(vValDonR.begin(), vValDonR.end());

  fprintf(stderr, "Probing SpliceMachine (starts)................");  
  fflush(stderr);

  strcpy(tempname,seqname);
  strcat(tempname,".spliceMSt");
  if (!ProbeFile(tempname)) SpliceMachine();

  ReadSMachineStarts(tempname, X->SeqLen);
  fprintf(stderr,"  done\n");

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
//  Read SpliceMachine Splice Site files
// -----------------------------
void SensorSMachine :: ReadSMachineSplices(char *name, int SeqLen)
{
  FILE *fp;
  char buf[FILENAME_MAX];
  char type[13];  // max len for "acceptor_rev"
  size_t len = 0;
  int i,pos,end =0;
  double force;

  if (!(fp = fopen(name, "r"))) {
    fprintf(stderr, "Cannot open splice sites file %s\n", name);
    exit(2);
  }
  // skip first line
  fgets(buf, FILENAME_MAX-1, fp);
  
  len =1;
  while (!end) {
    i = fscanf(fp,"%d %s %lf\n", &pos, type, &force);
    len ++;
    if (i == EOF) { end = 1; break;}
    if (i < 3) { end =2; break;}

    if (strcmp(type,"acceptor") == 0) {
      vPosAccF.push_back(pos);
      vValAccF.push_back(force);
    } else if (strcmp(type,"donor") == 0) {
      vPosDonF.push_back(pos-1);
      vValDonF.push_back(force);
    } else if (strcmp(type,"acceptor_rev") == 0) {
      vPosAccR.push_back(pos-1);
      vValAccR.push_back(force);
    } else if (strcmp(type,"donor_rev") == 0) {
      vPosDonR.push_back(pos);
      vValDonR.push_back(force);
    } else end = 2;
  }
  fclose(fp);

  if (end ==2) {
    fprintf(stderr, "Error in SpliceMachine splice site file %s, line %d\n", name, len);
    exit(2);
  }
}
// --------------------------
//  Read start file.
// --------------------------
void SensorSMachine :: ReadSMachineStarts(char *name, int Len)
{
  
  FILE *fp;
  char type[10];  // max len for "start_rev"
  int len = 0;
  int i,pos,end =0;
  double force;

  if (!(fp = fopen(name, "r"))) {
    fprintf(stderr, "Cannot open splice sites file %s\n", name);
    exit(2);
  }

  len =1;
  while (!end) {
    i = fscanf(fp,"%d %s %lf\n", &pos, type, &force);
    len ++;
    if (i == EOF) {end = 1; break;}
    if (i < 3) {end =2; break;}

    if (strcmp(type,"start") == 0) {
      vPosF.push_back(pos-1);
      vValF.push_back(force);
    } else if (strcmp(type,"start_rev") == 0) {
      vPosR.push_back(pos+1);
      vValR.push_back(force);
    } else end = 2;
  }
  fclose(fp);
  
  if (end ==2) {
    fprintf(stderr, "Error in SpliceMachine splice site file %s, line %d\n", name, len);
    exit(2);
  }
}

// ------------------------
//  GiveInfo signal SMachine. 
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
    PlotBarF(vPosF[i], (vPosF[i]%3)+1, 0.5, NORM(log(f),10.0), 2);
  }
  
  for (int i =0; i < (int)vPosR.size(); i++) {
    f = pow(vValR[i], startB)*(exp(-startP));
    PlotBarF(vPosR[i], -((X->SeqLen-vPosR[i])%3)-1, 0.5, NORM(log(f),10.0), 2);
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

  strcpy(s1,SMcmd);
  strcat(s1, " ");
  strcat(s1, PAR.getC("fstname"));

  system(s1);

  return;
}
