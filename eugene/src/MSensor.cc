#include "MSensor.h"

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
  for(int i=0; i<(int)theSensors.size(); i++) delete theSensors[i];
  theSensors.clear();

  for(int i=0; i<(int)msList.size(); i++) {
    delete msList[i];
    delete dllList[i];
  }
  msList.clear();
}

// ------------------------
//  Init Master.
// ------------------------
void MasterSensor :: InitMaster (void)
{
  char *key_name,*pluginsDir;
  int  val_prior;
  int i, j, nbSensors = 0;
  

  // On récupère les couples nom de sensor/priorité du .par
  PAR.ResetIter();
  while(PAR.getUseSensor(&key_name, &val_prior))
    msList.push_back(new UseSensor(val_prior, key_name));

  // On les tri
  sort(msList.begin(), msList.end(), Prior);
  
  // Liste des plugins devant se trouver dans ./PLUGINS
  soList  = new char*[(int)msList.size()];
  // Liste des plugins à utiliser (key pour attaquer PAR)
  useList = new char*[(int)msList.size()];
  pluginsDir = PAR.getC("EuGene.PluginsDir");

  for(i=0; i<(int)msList.size(); i++) {
    soList[i] = new char[FILENAME_MAX+1];
    strcpy(soList[i], pluginsDir);
    strcat(soList[i], msList[i]->Name);
    strcat(soList[i], ".so");
    useList[i] = new char[FILENAME_MAX+1];
    strcpy(useList[i], msList[i]->Name);
    strcat(useList[i], ".use");
  }
  
  for(i=0; i<(int)msList.size(); i++) 
    nbSensors += PAR.getI(useList[i]);
  
  dllList = new (SensorLoader *)[nbSensors];
  
  nbSensors = 0;

  for(i=0; i<(int)msList.size(); i++) 
    for (j=0; j<PAR.getI(useList[i]); j++) {   
      dllList[nbSensors] = new SensorLoader ( soList[i] );
      if(!dllList[nbSensors]->LastError()) {
	fprintf(stderr,"Loading %.21s",msList[i]->Name);
	for(int k=strlen(msList[i]->Name); k < 22; k++)
	  fprintf(stderr,".");
	fprintf(stderr,"%d\n",j);
	theSensors.push_back( dllList[nbSensors]->MakeSensor(j));
      }
      else {
	fprintf(stderr,"Error: ignored plugin (invalid or not found) : %s\n",soList[i]);
	exit(2);
      }
      nbSensors++;
    }

}

// ------------------------
//  Init. the sensors.
// ------------------------
void MasterSensor :: InitSensors (DNASeq *X)
{
  for(int i=0; i<(int)theSensors.size(); i++)
    theSensors[i]->Init(X);
}

// --------------------------
//  Get informations at pos.
// --------------------------
void MasterSensor :: GetInfoAt (DNASeq *X, int pos, DATA *d)
{
  int i;
  
  for(i=0; i <= DATA::Del;  i++)
    d->sig[i].Clear();
  
  for(i=0; i<= DATA::UTR3R; i++) d->contents[i] = 0.0;
  
  d->ESTMATCH_TMP = 0; // WARNING temporaire : EST -> on est dans intron
  
  for(i=0; i<(int)theSensors.size(); i++) 
    theSensors[i]->GiveInfo(X, pos, d);
  
  for (i=0; i<= DATA::Del; i++)
    d->sig[i].SetToDefault();
}

// --------------------------
//  Print informations at pos.
// --------------------------
void MasterSensor :: PrintDataAt (DNASeq *X, int pos, DATA *d)
{
  int i,j;
  printf("%6d %c", 1+pos, (*X)[pos]);

  for(i=0; i<=DATA::UTR3R; i++)
    printf (" %4.2f",d->contents[i]);

  for(i=0; i< Signal::LastEdge;  i++) {
    printf (" ||");
    for (j=0; j<= DATA::Del; j++)
      printf(" %4.2f",d->sig[j].weight[i]);
  }
  printf("\n");
}

// --------------------------------------------
//  Get special info at special pos.
//  Retourne TRUE si les sensors sont porteurs
//  d'infos de type "type" FALSE sinon.
// --------------------------------------------
int MasterSensor :: GetInfoSpAt (TYPE_SENSOR type,
				 DNASeq *X, int pos, DATA *d)
{
  int i;
  int info = FALSE;  // Aucune info
  
  for(i=0; i< DATA::Del+1;  i++)
    d->sig[i].Clear();
  
  for(i=0; i<DATA::UTR3R+1; i++) d->contents[i] = 0.0;
  
  d->ESTMATCH_TMP = 0; // WARNING temporaire : EST -> on est dans intron

  for(i=0; i<(int)theSensors.size(); i++) 
    if (theSensors[i]->type == type || theSensors[i]->type == Type_Multiple) {
      theSensors[i]->GiveInfo(X, pos, d);
      info = TRUE;
    }

  for (i=0; i<= DATA::Del; i++)
    d->sig[i].SetToDefault();

  return info;
}


// --------------------------
//  Post analyse the sensors.
// --------------------------
void MasterSensor :: PostAnalyse (Prediction *pred)
{
  for(int i=0; i<(int)theSensors.size(); i++)
    theSensors[i]->PostAnalyse(pred);
}

