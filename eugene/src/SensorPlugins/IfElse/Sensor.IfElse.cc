#include "Sensor.IfElse.h"
#include "../../EuGene/Dll.h"
#include "../../EuGene/MSensor.h"

extern Parameters PAR;
extern MasterSensor* MS;

/*************************************************************
 **                      SensorIfElse                       **
 *************************************************************/
// ----------------------
//  Default constructor.
// ----------------------
SensorIfElse :: SensorIfElse (int n, DNASeq *X) : Sensor(n)
{
  UseSensor *sensor_if, *sensor_else;
  int dll_index_if, dll_index_else;
  char *c = new char [FILENAME_MAX+1];

  type = Type_Multiple;

  // read sensors to use
  strcpy(c,"Sensor."); strcat(c, PAR.getC("IfElse.SensorIf",n));
  sensor_if = new UseSensor(0, c);
  strcpy(c,"Sensor."); strcat(c, PAR.getC("IfElse.SensorElse",n));
  sensor_else = new UseSensor(0, c);

  // load ".so" of sensors if necessary
  dll_index_if   = MS->LoadSensor(sensor_if,   n, "Sensor.IfElse :\n If   : ");
  dll_index_else = MS->LoadSensor(sensor_else, n, " Else : ");

  // create instance of sensors
  sensorIf   = MS->dllList[dll_index_if]->MakeSensor(n+1, X);
  sensorElse = MS->dllList[dll_index_else]->MakeSensor(n+1, X);
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
