#include "MSensor.h"
#include "SensorEuStop.h"
#include "SensorNStart.h"
#include "SensorNG2.h"
#include "SensorSPred.h"
#include "SensorMarkov.h"
#include "SensorBlastX.h"
#include "SensorRepeat.h"
#include "SensorRiken.h"
#include "SensorUser.h"
#include "SensorEst.h"



inline bool Prior(const UseSensor *A, const UseSensor *B)
{ return(A->Priority < B->Priority); }

/*************************************************************
 **                        UseSensor                        **
 *************************************************************/
// -------------------------
//  Default constructor.
// -------------------------
UseSensor :: UseSensor ()
{
  Priority = 0;
  Name[0] = '\000';
}

// -------------------------
//  Default constructor.
// -------------------------
UseSensor :: UseSensor (int priority, char name[FILENAME_MAX+1])
{
  Priority = priority;
  strcpy(Name,name);
}

// -------------------------
//  Default destructor.
// -------------------------
UseSensor :: ~UseSensor () {}



/*************************************************************
 **                       MasterSensor                      **
 *************************************************************/

extern Parameters PAR;

// ------------------------
// Default constructor.
// ------------------------
MasterSensor :: MasterSensor ()
{
}

// ------------------------
//  Default destructor.
// ------------------------
MasterSensor :: ~MasterSensor ()
{
  msList.clear();
  theSensors.clear();
}

// ------------------------
//  Init Master.
// ------------------------
void MasterSensor :: InitMaster ()
{
  char *key_name;
  int  val_prior;
  int i;
  
  while(PAR.getUseSensor(&key_name, &val_prior))
    msList.push_back(new UseSensor(val_prior, key_name));

  sort(msList.begin(), msList.end(), Prior);
  
  for(i=0;i<(int)msList.size();i++) {
    if(!strcmp(msList[i]->Name, "sensor.EuStop") && msList[i]->Priority != 0) {
      theSensors.push_back(new SensorEuStop());
    }
    else if(!strcmp(msList[i]->Name, "sensor.NStart") && msList[i]->Priority != 0) {
      theSensors.push_back(new SensorNStart());
    }
    else if(!strcmp(msList[i]->Name, "sensor.NG2") && msList[i]->Priority != 0) {
      theSensors.push_back(new SensorNG2());
    }
    else if(!strcmp(msList[i]->Name, "sensor.SPred") && msList[i]->Priority != 0) {
      theSensors.push_back(new SensorSPred());
    }
    else if(!strcmp(msList[i]->Name, "sensor.Markov") && msList[i]->Priority != 0) {
      theSensors.push_back(new SensorMarkov());
    }
    else if(!strcmp(msList[i]->Name, "sensor.BlastX") && msList[i]->Priority != 0) {
      if(PAR.getI("blastopt"))
	theSensors.push_back(new SensorBlastX());
    }
    else if(!strcmp(msList[i]->Name, "sensor.Repeat") && msList[i]->Priority != 0) {
      if(PAR.getI("ncopt"))
	theSensors.push_back(new SensorRepeat());
    }
    else if(!strcmp(msList[i]->Name, "sensor.Riken") && msList[i]->Priority != 0) {
      if(PAR.getI("raflopt"))
	theSensors.push_back(new SensorRiken());
    }
    else if(!strcmp(msList[i]->Name, "sensor.User") && msList[i]->Priority != 0) {
      if(PAR.getI("userinfo"))
	theSensors.push_back(new SensorUser());
    }
    else if(!strcmp(msList[i]->Name, "sensor.Est") && msList[i]->Priority != 0) {
      if(PAR.getI("estopt"))
	theSensors.push_back(new SensorEst());
    }
    else if(msList[i]->Priority == 0)
      fprintf(stderr,"WARNING: Ignored information : %s\n", msList[i]->Name);
    else
      fprintf(stderr,"WARNING: Unknown information : %s\n", msList[i]->Name);
  }
  msList.clear();
}

// ------------------------
//  Init. the sensors.
// ------------------------
void MasterSensor :: InitSensors (DNASeq *X)
{
  for(int i=0; i<(int)theSensors.size(); i++) {
    theSensors[i]->Init(X);
  }
}

// --------------------------
//  Get informations at pos.
// --------------------------
void MasterSensor :: GetInfoAt (DNASeq *X, int pos, DATA *d)
{
  int i;
  for(i=0; i<2;  i++) d->Stop[i] = d->Start[i] = d->Acc[i] = d->Don[i] = 0.0;
  for(i=0; i<13; i++) d->ContentScore[i] = 0.0;

  d->ESTMATCH_TMP = 0; // WARNING temporaire : EST -> on est dans intron

  for(i=0; i<(int)theSensors.size(); i++) {
    theSensors[i]->GiveInfo(X, pos, d);
  }
}

// --------------------------------------------
//  Get special info at pos.
//  Retourne TRUE si les sensors sont porteurs
//  d'infos de type "type" FALSE sinon.
// --------------------------------------------
int MasterSensor :: GetInfoSpAt (TYPE_SENSOR type,
				 DNASeq *X, int pos, DATA *d)
{
  int i;
  int info = FALSE;  // Aucune info

  for(i=0; i<2;  i++) d->Stop[i] = d->Start[i] = d->Acc[i] = d->Don[i] = 0.0;
  for(i=0; i<13; i++) d->ContentScore[i] = 0.0;

  d->ESTMATCH_TMP = 0; // WARNING temporaire : EST -> on est dans intron

  for(i=0; i<(int)theSensors.size(); i++) {
    if(theSensors[i]->type == type || theSensors[i]->type == Type_Multiple) {
      theSensors[i]->GiveInfo(X, pos, d);
      info = TRUE;
    }
    else if(theSensors[i]->type == Type_Unknown)
      return info;  // Aucune info pour ce type
  }
  return info;
}

// --------------------------
//  Reset type.
// --------------------------
void MasterSensor :: ResetType ()
{
  for(int i=0; i<(int)theSensors.size(); i++)
    theSensors[i]->type = Type_Unknown;
}

// -------------------------
//  Destroy the sensors.
// -------------------------
void MasterSensor :: ResetSensors ()
{
  for(int i=0; i<(int)theSensors.size(); i++)
    delete theSensors[i];
}
