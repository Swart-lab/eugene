#include "MSensor.h"


// reserve memory for static variables
bool MasterSensor::IsInitialized = false;
std::string PluginsDir;
std::vector <UseSensor*> MasterSensor::LoadedSensorsList;
std::vector <std::string> MasterSensor::soList;
std::vector <std::string> MasterSensor::useList;
std::vector <SensorLoader*> MasterSensor::dllList;


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
//  Default destructor.
// ------------------------
MasterSensor :: ~MasterSensor (void)
{
  ///////////////////  Need to comment ///////////////////////
  // Create a segmentation fault when asking to delete a second sensor instance 
  // of sensor with static attributs; for example MarkovIMM.
  //  for(int i=0; i<(int)theSensors.size(); i++) delete theSensors[i];
  //  theSensors.clear();
  ////////////////////////////////////////////////////////////
  //  for(int i=0; i<(int)msList.size(); i++) delete msList[i];
  //  msList.clear();
}

// ------------------------
//  Init Master.
// ------------------------
void MasterSensor :: InitMaster (DNASeq *X)
{
  char *key_name;
  int  val_prior;
  int i, j;
  char* c = new char[FILENAME_MAX+1];
  vector <UseSensor*> msList;

  if (!IsInitialized) {
    PluginsDir = PAR.getC("EuGene.PluginsDir");

    // On récupère les couples nom de sensor/priorité du .par
    PAR.ResetIter();
    while(PAR.getUseSensor(&key_name, &val_prior)) 
      msList.push_back(new UseSensor(val_prior, key_name));

    // On les tri
    sort(msList.begin(), msList.end(), Prior);

    // Load used sensors
    for(i=0; i<(int)msList.size(); i++) {
      strcpy(c, msList[i]->Name); 
      strcat(c, ".use");
      for (j=0; j<PAR.getI(c); j++)
	LoadSensor( msList[i], j, "" );
    }

    IsInitialized = true;
  }

  // create instance(s) of used sensors
  for(i=0; i<(int)LoadedSensorsList.size(); i++) {
    strcpy(c, useList[i].c_str()); // to avoid passing a const as argument to PAR.getI
    for (j=0; j<PAR.getI(c); j++) 
      theSensors.push_back( dllList[i]->MakeSensor(j, X));
  }
}

// ------------------------
// load the sensor s  if no yet done (j is just for print)
// return index of sensor in dllList
// ------------------------
int MasterSensor :: LoadSensor (UseSensor* s, int j, char *stdErr)
{
  bool is_loaded = false;
  int dll_index = 0;

  // is the sensor yet loaded
  for (unsigned int i=0; i< LoadedSensorsList.size(); i++)
    if ( ((string) s->Name) == ((string) LoadedSensorsList[i]->Name) ) {
      is_loaded = true;
      dll_index = i;
      i = LoadedSensorsList.size();
    }

  if (!is_loaded) {
    fprintf(stderr, "%s", stdErr);
    soList.push_back ( PluginsDir + s->Name + ".so" );
    useList.push_back ( ((string) s->Name) + ".use" );
    LoadedSensorsList.push_back( s );

    dll_index = soList.size() - 1; 

    dllList.push_back( new SensorLoader ( soList[dll_index].c_str() ) );
    if(!dllList[dll_index]->LastError()) {
      fprintf(stderr,"Loading %.21s", s->Name);
      for(int k=strlen(s->Name); k < 22; k++) fprintf(stderr,".");
      fprintf(stderr,"%d\n",j);
    } else {
      fprintf(stderr,"Error: ignored plugin (invalid or not found) : %s\n",soList[dll_index].c_str());
      exit(2);
    }
  }
  return dll_index;
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

