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
// File:     Sensor.NG2.cc
// Contents: Sensor NG2
// ------------------------------------------------------------------

#include "Sensor.NG2.h"

extern Parameters PAR;

#define NORM(x,n) (((n)+(Max(-(n),x)))/(n))

/*************************************************************
 **                        SensorNetGene2                   **
 *************************************************************/

// ----------------------
//  Default constructor.
// ----------------------
SensorNG2 :: SensorNG2 (int n, DNASeq *X) : Sensor(n)
{
  char tempname[FILENAME_MAX+1];

  type = Type_Acc|Type_Don;
    
  fprintf(stderr, "Reading splice site file (NetGene2)...........");  
  fflush(stderr);
  strcpy(tempname,PAR.getC("fstname"));
  strcat(tempname,".splices");
  ReadNG2F(tempname, X->SeqLen);
  fprintf(stderr,"forward,");
  fflush(stderr);

  strcpy(tempname,PAR.getC("fstname"));
  strcat(tempname,".splicesR");
  ReadNG2R(tempname, X->SeqLen);
  fprintf(stderr," reverse done\n");
  
  CheckSplices(X,vPosAccF, vPosDonF, vPosAccR, vPosDonR);
  
  // vectors for reverse are put in the increasing order
  reverse(vPosAccR.begin(), vPosAccR.end()); 
  reverse(vValAccR.begin(), vValAccR.end());
  reverse(vPosDonR.begin(), vPosDonR.end());
  reverse(vValDonR.begin(), vValDonR.end());
}

// ----------------------
//  Default destructor.
// ----------------------
SensorNG2 :: ~SensorNG2 ()
{
  vPosAccF.clear();  vPosAccR.clear();
  vPosDonF.clear();  vPosDonR.clear();
  vValAccF.clear();  vValAccR.clear();
  vValDonF.clear();  vValDonR.clear();
}

// ----------------------
//  Init NG2.
// ----------------------
void SensorNG2 :: Init (DNASeq *X)
{
  accB = PAR.getD("NG2.accB*",GetNumber());
  accP = PAR.getD("NG2.accP*",GetNumber());
  donB = PAR.getD("NG2.donB*",GetNumber());
  donP = PAR.getD("NG2.donP*",GetNumber());

  iAccF = iDonF = iAccR = iDonR = 0;
  PositionGiveInfo = -1;

  if (PAR.getI("Output.graph")) Plot(X);
}

// -----------------------------
//  Read NetGene2 forward file.
// -----------------------------
void SensorNG2 :: ReadNG2F(char name[FILENAME_MAX+1], int SeqLen)
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
void SensorNG2 :: ReadNG2R(char name[FILENAME_MAX+1], int SeqLen)
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

// ------------------------
//  GiveInfo signal NG2.   pow(atof(sdon),  donB)*donP
// ------------------------
void SensorNG2 :: GiveInfo (DNASeq *X, int pos, DATA *d)
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

}

// ----------------------------
//  Plot Sensor information
// ----------------------------
void SensorNG2 :: Plot(DNASeq *X)
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
}

// ------------------
//  Post analyse
// ------------------
void SensorNG2 :: PostAnalyse(Prediction *pred)
{
}
