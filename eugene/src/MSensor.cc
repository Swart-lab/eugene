#include "MSensor.h"


// reserve memory for static variables
bool                        MasterSensor::IsInitialized = false;
std::string                 MasterSensor::PluginsDir;
std::vector <std::string>   MasterSensor::LoadedSensorsList;
std::vector <std::string>   MasterSensor::MSSensorsList;
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




/*************************************************************
 **                       MasterSensor                      **
 *************************************************************/

extern Parameters PAR;


// ------------------------
//  Default destructor.
// ------------------------
MasterSensor :: ~MasterSensor (void)
{
  // Need to be commented because a segmentation fault sometimes occured 
  // for example in Units tests
  // delete created instances with new
  //unsigned int i;
  //for(i=0; i<theSensors.size(); i++) delete theSensors[i];
  //theSensors.clear();
  //for(i=0; i<dllList.size(); i++) delete dllList[i];
  //dllList.clear();
}

// ------------------------
//  Init Master.
// ------------------------
void MasterSensor :: InitMaster (DNASeq *X)
{
  char *key_name;
  int  val_prior;
  int  j;
  unsigned int i;
  char* c = new char[FILENAME_MAX+1];
  std::vector <UseSensor*> msList;
  std::string use_name;

  if (!IsInitialized) {
    PluginsDir = PAR.getC("EuGene.PluginsDir");

    // On récupère les couples nom de sensor/priorité du .par
    PAR.ResetIter();
    while( PAR.getUseSensor(&key_name, &val_prior) ) 
      msList.push_back( new UseSensor(val_prior, key_name) );

    // On les tri
    sort(msList.begin(), msList.end(), Prior);
    
    // Update the list of sensors defined at the top level
    for (i=0; i<msList.size(); i++)  MSSensorsList.push_back( (std::string) msList[i]->Name );

    // delete instances of UseSensor
    for (i=0; i<msList.size(); i++)  delete msList[i];
    
    IsInitialized = true;
  } 

  // create instance(s) of used sensors
  for (i=0; i<MSSensorsList.size(); i++) {
    use_name = MSSensorsList[i] + ".use";
    strcpy(c, use_name.c_str()); // to avoid passing a const as argument to PAR.getI
    for (j=0; j<PAR.getI(c); j++) 
      theSensors.push_back( MakeSensor( MSSensorsList[i].c_str(),j, X) );
  }
  delete [] c;
}


// ------------------------
// Create an instance of a sensor and load the .so before if necessary
// ------------------------
Sensor* MasterSensor :: MakeSensor (std::string name, int n, DNASeq *X)
{
  Sensor* s = NULL;
  int dll_index;

  dll_index = LoadSensor (name);
  s = dllList[dll_index]->MakeSensor( n, X);

  return s;
}


// ------------------------
// load the sensor s  if no yet done (j is just for print)
// return index of sensor in dllList
// ------------------------
int MasterSensor :: LoadSensor (std::string name)
{
  bool is_loaded = false;
  int dll_index = 0;
  std::string complete_name;

  // is the sensor yet loaded
  for (unsigned int i=0; i< LoadedSensorsList.size(); i++)
    if ( name == LoadedSensorsList[i] ) {
      is_loaded = true;
      dll_index = i;
      i = LoadedSensorsList.size();
    }

  if (!is_loaded) {
    LoadedSensorsList.push_back( name );
    dll_index = LoadedSensorsList.size() - 1; 
    complete_name = PluginsDir + name + ".so";
    dllList.push_back( new SensorLoader ( complete_name.c_str() ) );
    if(!dllList[dll_index]->LastError()) {
      fprintf(stderr,"Loading %.21s", name.c_str());
      for(int k=name.size(); k < 22; k++) fprintf(stderr,".");
      fprintf(stderr,"done\n");
    } else {
      fprintf(stderr,"Error: ignored plugin (invalid or not found) : %s\n",complete_name.c_str());
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
int MasterSensor :: GetInfoSpAt (unsigned char type,
				 DNASeq *X, int pos, DATA *d)
{
  int i;
  int info = FALSE;  // Aucune info
  
  for(i=0; i< DATA::Del+1;  i++)
    d->sig[i].Clear();
  
  for(i=0; i<DATA::UTR3R+1; i++) d->contents[i] = 0.0;
  
  d->ESTMATCH_TMP = 0; // WARNING temporaire : EST -> on est dans intron

  for(i=0; i<(int)theSensors.size(); i++) 
    if (theSensors[i]->type & type) {
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

