/*****************************************************************************/
/*             Copyright (c) 2004 by INRA. All rights reserved.              */
/*                 Redistribution is not permitted without                   */
/*                 the express written permission of INRA.                   */
/*                     Mail : tschiex@toulouse.inra.fr                       */
/*---------------------------------------------------------------------------*/
/* File         : EuGeneTk/SensorPlugins/PatConst/Sensor.PatConst.cc         */
/* Description  : Sensor pattern const                                       */
/* Authors      : P.Bardou, S.Foissac, M.J.Cros, A.Moisan, T.Schiex          */
/* History      : Fev 2004         	   		                     */
/*****************************************************************************/

#include "Sensor.PatConst.h"

/*************************************************************
 **                       SensorPatConst
 *************************************************************/

extern Parameters PAR;

// ----------------------
// Default constructor.
// ----------------------
SensorPatConst :: SensorPatConst (int n) : Sensor(n)
{
  newStatePos = PAR.getI("PatConst.newStatePos", GetNumber());
  patType     = PAR.getC("PatConst.type",        GetNumber());
  pattern     = PAR.getC("PatConst.pat",         GetNumber());
  patLen      = strlen(pattern);
  for(int i=0; i<patLen; i++)
    pattern[i] = tolower(pattern[i]);
  
  // Type initialisation

  if(!strcmp(patType, "start")) {
    sigTypeIndex = DATA::Start;   type = Type_Start;  }
  else if(!strcmp(patType, "insertion")) {
    sigTypeIndex = DATA::Ins;     type = Type_FS; }
  else if(!strcmp(patType, "deletion"))  {
    sigTypeIndex = DATA::Del;     type = Type_FS; }
  else if(!strcmp(patType, "transstart")) {
    sigTypeIndex = DATA::tStart;  type = Type_TStart; }
  else if(!strcmp(patType, "transstop")) {
    sigTypeIndex = DATA::tStop;   type = Type_TStop; }
  else if(!strcmp(patType, "donor")) {
    sigTypeIndex = DATA::Don;     type = Type_Don; }
  else if(!strcmp(patType, "acceptor")) {
    sigTypeIndex = DATA::Acc;     type = Type_Acc; }
  else if(!strcmp(patType, "stop")) {
    sigTypeIndex = DATA::Stop;     type = Type_Stop;}
  else {
    fprintf(stderr, "Error \"%s\" undefined type in the parameter file"
	    " for PatConst.pat[%d]\n", patType, GetNumber());
    fprintf(stderr, "Type must be : start, insertion, deletion, "
	    "transstart, transstop, stop,\n               "
	    "acceptor or donor.\n");
    exit(2);
  }
}

// ----------------------
//  Default destructor.
// ----------------------
SensorPatConst :: ~SensorPatConst ()
{
}
// ----------------------
//  Init.
// ----------------------
void SensorPatConst :: Init (DNASeq *X)
{
  patP   = PAR.getD("PatConst.patP*",  GetNumber());
  patPNo = PAR.getD("PatConst.patPNo*",GetNumber());

  //if (PAR.getI("Output.graph")) Plot(X);
}

// ------------------------
//  GiveInfo.
// ------------------------
void SensorPatConst :: GiveInfo (DNASeq *X, int pos, DATA *d)
{
  int  i;
  bool isSigF = true;
  bool isSigR = true;
 
  for(i=0; i<patLen; i++)
    if((*X)[pos+i-newStatePos+1] != pattern[i]) {
      isSigF=false;
      break;
    }
  for(i=patLen-1; i>-1; i--)
    if((*X)(pos-i+newStatePos-2) != pattern[i]) {
      isSigR=false;
      break;
    }
  
  if(isSigF) {
    d->sig[sigTypeIndex].weight[Signal::Forward]   += patP;
    d->sig[sigTypeIndex].weight[Signal::ForwardNo] += patPNo;
  }
  if(isSigR) {
    d->sig[sigTypeIndex].weight[Signal::Reverse]   += patP;
    d->sig[sigTypeIndex].weight[Signal::ReverseNo] += patPNo;
  }
}

// ----------------------------
//  Plot Sensor information.
// ----------------------------
void SensorPatConst :: Plot(DNASeq *X)
{
}

// ------------------
//  Post analyse.
// ------------------
void SensorPatConst :: PostAnalyse(Prediction *pred)
{
}
