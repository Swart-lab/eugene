#include "Sensor.IfElse.h"
#include "../../EuGene/Dll.h"

extern Parameters PAR;

/*************************************************************
 **                      SensorIfElse                       **
 *************************************************************/
// ----------------------
//  Default constructor.
// ----------------------
SensorIfElse :: SensorIfElse (int n, DNASeq *X) : Sensor(n)
{
  char * sensorName = new char[FILENAME_MAX+1];
  char * pluginsDir;

  type = Type_Multiple;

  pluginsDir = PAR.getC("EuGene.PluginsDir");

  strcpy(sensorName, pluginsDir);
  strcat(sensorName, "Sensor.");
  strcat(sensorName, PAR.getC("IfElse.SensorIf",n));
  strcat(sensorName, ".so");
  sensorLIf =  new SensorLoader (sensorName);
  if (!sensorLIf->LastError()) {
    fprintf(stderr," -If  : Sensor.%s\t%d\n",PAR.getC("IfElse.SensorIf",n),n+1);
    sensorIf = sensorLIf->MakeSensor(n+1, X);
  }
  else {
    fprintf(stderr,"WARNING: ignored plugin (invalid or not found) : %s\n",
	    PAR.getC("IfElse.SensorIf"));
    exit(2);
  }
  
  strcpy(sensorName, pluginsDir);
  strcat(sensorName, "Sensor.");
  strcat(sensorName, PAR.getC("IfElse.SensorElse",n) );
  strcat(sensorName, ".so");
  sensorLElse = new SensorLoader (sensorName);
  if (!sensorLElse->LastError()) {
    fprintf(stderr," -Else: Sensor.%s\t%d\n",PAR.getC("IfElse.SensorElse",n),n+1);
    sensorElse = sensorLElse->MakeSensor(n+1, X);
  }
  else { 
    fprintf(stderr,"WARNING: ignored plugin (invalid or not found) : %s\n",
	    PAR.getC("IfElse.SensorElse"));
    exit(2);
  }
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

//
void MergeInfo(DATA *ifData, DATA *elseData, DATA *d)
{
  int i,j;

  // on fusionne les signaux dans ifData
  for (j=0; j < DATA::LastSigType; j++) {
    for (i=0; i<= Signal::Reverse;  i++)
      if (!ifData->sig[j].IsSet(i)) {
	ifData->sig[j].weight[i] = elseData->sig[j].weight[i];
	ifData->sig[j].weight[i+2] = elseData->sig[j].weight[i+2];
      }
  }

  // on repercute dans d
  for (j=0; j < DATA::LastSigType; j++) {
    for (i=0; i<= Signal::Reverse;  i++)
      if (ifData->sig[j].IsSet(i)) {
	d->sig[j].weight[i] += ifData->sig[j].weight[i];
	d->sig[j].weight[i+2] += ifData->sig[j].weight[i+2];
      }
  }

  for(i=0; i< DATA::LastContentsType; i++)
    d->contents[i] +=
      (ifData->contents[i] ? ifData->contents[i] : elseData->contents[i]);
}

// -----------------------
//  GiveInfo signal stop.
// -----------------------
void SensorIfElse :: GiveInfo (DNASeq *X, int pos, DATA *d)
{
  DATA ifData;
  DATA elseData;
  int i;

  for(i=0; i<  DATA::LastSigType;  i++)
    ifData.sig[i].Clear();
  for(i=0; i< DATA::LastContentsType; i++) ifData.contents[i] = 0.0;

  sensorIf->GiveInfo(X,pos,&ifData);

  for(i=0; i<  DATA::LastSigType;  i++)
    elseData.sig[i].Clear();
  for(i=0; i< DATA::LastContentsType; i++) elseData.contents[i] = 0.0;


  sensorElse->GiveInfo(X,pos,&elseData);
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

// ------------------
//  Post analyse
// ------------------
void SensorIfElse :: PostAnalyse(Prediction *pred)
{
  sensorIf->PostAnalyse(pred);
  sensorElse->PostAnalyse(pred);
}
