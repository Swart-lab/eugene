#include "Sensor.IfElse.h"
#include "../Dll.h"

extern Parameters PAR;

/*************************************************************
 **                      SensorIfElse                       **
 *************************************************************/
// ----------------------
//  Default constructor.
// ----------------------
SensorIfElse :: SensorIfElse (int n) : Sensor(n)
{
  char * sensorName;

  sensorName = new char[FILENAME_MAX+1];
  strcpy(sensorName, "./PLUGINS/Sensor.");
  strcat(sensorName, PAR.getC("IfElse.SensorIf"));
  strcat(sensorName, ".so");
  sensorLIf =  new SensorLoader (sensorName);

  if (!sensorLIf->LastError()) {
    fprintf(stderr,"Loading Sensor.%s, %d\n",PAR.getC("IfElse.SensorIf"),0);
    sensorIf = sensorLIf->MakeSensor(0);
  } else 
    fprintf(stderr,"WARNING: ignored plugin (invalid or not found) : %s\n",
	    PAR.getC("IfElse.SensorIf"));

  strcpy(sensorName, "./PLUGINS/Sensor.");
  strcat(sensorName, PAR.getC("IfElse.SensorElse"));
  strcat(sensorName, ".so");
  sensorLElse =  new SensorLoader (sensorName);
  if (!sensorLElse->LastError()) {
    fprintf(stderr,"Loading Sensor.%s, %d\n",PAR.getC("IfElse.SensorElse"),0);
    sensorElse = sensorLElse->MakeSensor(0);
  } else 
    fprintf(stderr,"WARNING: ignored plugin (invalid or not found) : %s\n",
	    PAR.getC("IfElse.SensorElse"));
}
// ----------------------
//  Default destructor.
// ----------------------
SensorIfElse :: ~SensorIfElse ()
{
  delete sensorLIf;
  delete sensorLElse;
}
// ----------------------
//  Init from new seq.
// ----------------------
void SensorIfElse :: Init (DNASeq *X)
{
  sensorIf->Init(X);
  sensorElse->Init(X);
}

// -----------------------
//  ResetIter.
// -----------------------
void SensorIfElse :: ResetIter ()
{
  sensorIf->ResetIter();
  sensorElse->ResetIter();
}

//
void MergeInfo(DATA *ifData, DATA *elseData, DATA *d)
{
  int i;

  // on fusionne les signaux dans ifData
  for(i=0; i<2;  i++) {
    if (!ifData->Stop[i]) ifData->Stop[i] = elseData->Stop[i];
    if (!ifData->Start[i]) ifData->Start[i] = elseData->Start[i];
    if (!ifData->Acc[i]) ifData->Acc[i] = elseData->Acc[i];
    if (!ifData->Don[i]) ifData->Don[i] = elseData->Don[i];
  }
  
// on repercute dans d
  for(i=0; i<2;  i++) {
    if (ifData->Stop[i]) d->Stop[i] = ifData->Stop[i];
    if (ifData->Start[i]) d->Start[i] = ifData->Start[i];
    if (ifData->Acc[i]) d->Acc[i] = ifData->Acc[i];
    if (ifData->Don[i]) d->Don[i] = ifData->Don[i];
  }

 for(i=0; i<13; i++) 
   d->ContentScore[i] += 
     (ifData->ContentScore[i] ? ifData->ContentScore[i] : elseData->ContentScore[i]);    
}

// -----------------------
//  GiveInfo signal stop.
// -----------------------
void SensorIfElse :: GiveInfo (DNASeq *X, int pos, DATA *d)
{
  DATA ifData;
  DATA elseData;
  int i;

  for(i=0; i<2;  i++) 
    ifData.Stop[i] = ifData.Start[i] = ifData.Acc[i] = ifData.Don[i] = 0.0;
  for(i=0; i<13; i++) 
    ifData.ContentScore[i] = 0.0;

  sensorIf->GiveInfo(X,pos,&ifData);

  for(i=0; i<2;  i++) 
    elseData.Stop[i] = elseData.Start[i] = elseData.Acc[i] = elseData.Don[i] = 0.0;
  for(i=0; i<13; i++) 
    elseData.ContentScore[i] = 0.0;

  sensorElse->GiveInfo(X,pos,&elseData);
  MergeInfo(&ifData,&elseData,d);
}
// -------------------------
//  GiveInfoAt signal stop.
// -------------------------
void SensorIfElse :: GiveInfoAt (DNASeq *X, int pos, DATA *d)
{
  DATA ifData;
  DATA elseData;
  int i;

  for(i=0; i<2;  i++) 
    ifData.Stop[i] = ifData.Start[i] = ifData.Acc[i] = ifData.Don[i] = 0.0;
  for(i=0; i<13; i++) 
    ifData.ContentScore[i] = 0.0;

  sensorIf->GiveInfoAt(X,pos,&ifData);

  for(i=0; i<2;  i++) 
    elseData.Stop[i] = elseData.Start[i] = elseData.Acc[i] = elseData.Don[i] = 0.0;
  for(i=0; i<13; i++) 
    elseData.ContentScore[i] = 0.0;

  sensorElse->GiveInfoAt(X,pos,&elseData);
  MergeInfo(&ifData,&elseData,d);
}
// ----------------------------
//  Plot Sensor information
// ----------------------------
void SensorIfElse :: Plot(DNASeq *X)
{
  sensorIf->Plot(X);
  sensorElse->Plot(X);
}
