#include "Param.h"

// -----------------------------------------------
//  Gestion des arguments (integer)
// -----------------------------------------------
int TestIArg(char *arg)
{
  int tmp;
  return (sscanf(arg, "%d", &tmp) == 1) && (tmp >= 0);
}

// -----------------------------------------------
//  Gestion des arguments (double)
// -----------------------------------------------
int TestDArg(char *arg)
{
    double tmp;
    return  (sscanf(arg, "%lf", &tmp) == 1) && (tmp >= 0);
}

// ------------------------
//  Default constructor.
// ------------------------
Parameters :: Parameters ()
{
}

// ------------------------
//  initParam.
// ------------------------
void Parameters :: initParam (int argc, char * argv[])
{
  fprintf(stderr,"EuGene rel. %s\n",VERSION);
  fprintf(stderr,"Loading parameters file...");
  fflush(stderr);

  ReadPar(argv[0]);
  ReadArg(argc, argv);

  if (m.count("Param.debug")) 
    for (iter = m.begin(); iter!=m.end(); ++iter)
      fprintf(stderr,"%s = %s\n",iter->first, iter->second);

  iter = m.begin();  // Cf. : getUseSensor
}

// ------------------------
//  Read arguments line.
// ------------------------
void Parameters :: ReadArg(int argc, char * argv[])
{
  int carg, errflag = 0;
  char *key, *val = NULL;
  char* indexPos = NULL;
  
  m["fstname"] = "\000";        // no default input

  while ((carg = getopt(argc, argv, "UREdrshm:w:f:n:o:p:x:y:c:u:v:g::b::l:O:D:t::")) != -1) {
    
    switch (carg) {
      
    case 'D':           /* Definition of any parameter */

      if ( (indexPos = index(optarg,'=')) && (indexPos > optarg) && (indexPos < optarg+strlen(optarg)-1)) {
	key = new char[FILENAME_MAX+1];
	val = new char[FILENAME_MAX+1];
	*indexPos = 0;
	sscanf(optarg, "%s", key);
	sscanf(indexPos+1, "%s", val);
	m[key] = val;
      }
      else {
	fprintf(stderr,"Invalid command line parameter definition: %s\n",optarg);
	exit(2);
      }
      break;

    case 'n':           /* -n normalize across frames      */
      if (!TestIArg(optarg))
	errflag++;
      else m["Output.normopt"] = optarg;
      break;
      
    case 'h':           /* help                             */
      errflag++;
      break;
      
    case 's':           /* Single gene mode. Start/End on IG only   */
      m["EuGene.ExonPrior"]   = "0.0";
      m["EuGene.IntronPrior"] = "0.0";
      m["EuGene.FivePrimePrior"]  = "0.0";
      m["EuGene.ThreePrimePrior"] = "0.0";
      break;

    case 'r':    /* RepeatMasker input */
      m["Sensor.Repeat.use"] = "TRUE";
      break;
      
    case 'c':           /* -c couverture      */
      if (! TestIArg(optarg))
        errflag++;
      else m["Output.golap"] = optarg;
      break;
      
    case 'l':           /* -l imglen      */
      if (! TestIArg(optarg))
        errflag++;
      else m["Output.glen"] = optarg;
      break;
      
    case 'u':           /* -u From      */
      if (! TestIArg(optarg))
	errflag++;
      else m["Output.gfrom"] = optarg;
      break;
      
    case 'v':           /* -v To      */
      if (! TestIArg(optarg))
	errflag++;
      else m["Output.gto"] = optarg;
      break;
      
    case 'x':           /* -x resx      */
      if (! TestIArg(optarg))
        errflag++;
      else m["Output.resx"] = optarg;
      break;
      
    case 'y':           /* -y resy      */
      if (! TestIArg(optarg))
        errflag++;
      else m["Output.resy"] = optarg;
      break;
      
    case 'g':          /* -g "Graphic File"                */
      if (optarg) {
	if (strspn(".fasta", optarg)    != 6 
	    && strspn(".fsa", optarg)   != 4
	    && strspn(".tfa", optarg)   != 4
	    && strspn(".txt", optarg)   != 4
	    && optarg[0] != '-') {
	  if(strcmp(optarg, "NO"))     // Argument suivant le -g
	    m["grnameArg"] = optarg;   // != NO pris en compte
	  else 
	    m["grnameArg"] = "0";      // == NO non pris en compte
	  m["Output.graph"] = "TRUE";
	}
	else errflag++;
      }
      else {
	m["Output.graph"] = "TRUE";
	m["grnameArg"] = "0";
      }
      break;
      
    case 'p':           /* print opt: short/long/detailed   */
      m["Output.format"] = optarg;
      if ((optarg[0] == 'h') &&  getI("Output.graph") == 0) {
	m["Output.graph"] = "TRUE"; // HTML output means graphical output 
	m["grnameArg"]    = "0";
      }
      if ((optarg[0] != 's') && (optarg[0] != 'l') && (optarg[0] != 'g') &&
	  (optarg[0] != 'd') && (optarg[0] != 'h') && (optarg[0] != 'a'))
	errflag++;
      break;
      
    case 'm':           /* -m matrix                        */
      m["Markov.matname"] = optarg;
      break;

    case 'O':           /* -O output                        */
      m["Output.Prefix"] = optarg;
      break;

    case 'o':           /* -o offset                        */
      if (! TestIArg(optarg))
	errflag++;
      else m["Output.offset"] = optarg;
      break;
      
    case 'w':           /* -w window                        */
      if (! TestIArg(optarg))
	errflag++;
      else m["Output.window"] = optarg;
      break;

    case 'b':           /* -b  use blastx result            */
      m["Sensor.BlastX.use"] = "TRUE";
      if (optarg) {
	m["BlastX.levels"] = optarg;
	if(strlen(optarg) > 8)
	  errflag++;
      }
      else m["BlastX.levels"] = "123";
      break;
      
    case 'E':
      m["Est.PostProcess"] = "TRUE";
      break;
      
    case 'U':
      m["Sensor.User.use"] = "TRUE";
      break;

    case 'd':           /* -d  use cDNA blastn results      */
      m["Sensor.Est.use"] = "TRUE";
      break; 

    case 'R':           /* -R use RAFL-like EST*/
      m["Sensor.Riken.use"] = "TRUE";
      break;

    case 'f':           /* -f frameshift loglike            */
      if (! TestDArg(optarg))
	errflag++;
      else m["EuGene.frameshift"] = optarg;

      break;
      
    case 't':           // -t use tblastx results   
      m["Sensor.Homology.use"] = "TRUE";
      if (optarg) {
	m["Homology.protmatname"] = optarg;
      }
      break;

    case '?':           /* bad option                       */
      errflag++;
    }
  }

  // may remain arguments -> fasta filenames
  if ((argc - optind) == 0)
    errflag++;
  
  // check usage
  if (errflag) {
    fprintf(stderr, "Usage: EuGene [-h] [-m matrix] [-n 0|1|2] [-s] [-p h|g|s|l|d|a]\n"
	    "              [-w window] [-b {levels}] [-d] [-R] [-E] [-U] [-o offset]\n"
	    "              [-g {graphArg}] [-u start] [-v end] [-l len] [-c olap]\n"
	    "              [-x xres] [-y yres] FASTA files\n");
    exit(1);
  }
}

// ------------------------
//  Read parameters file.
// ------------------------ 
void Parameters :: ReadPar(char *argv)
{
  char line    [MAX_LINE];
  char tempname[FILENAME_MAX+1];
  char *key, *val = NULL;
  
  strcpy(tempname, argv);
  strcat(tempname, ".par");
  fp = FileOpen(EugDir, BaseName(tempname),"r");

  while(fgets (line, MAX_LINE, fp) != NULL) {
    if (line[0] != '#') {
      key = new char[FILENAME_MAX+1];
      val = new char[FILENAME_MAX+1];
      if(sscanf(line, "%s %s", key, val) == 2)  
	m[key] = val;
      else {
	if (EugDir != NULL)
	  fprintf(stderr, "Incorrect parameter file %s/%s\n",EugDir,BaseName(tempname));
	else
	  fprintf(stderr, "Incorrect parameter file %s\n",BaseName(tempname));
	exit(2);
      }
    }
  }

  fclose(fp);
  fprintf(stderr,"done ");
  if(strcmp(getC("EuGene.versionPAR"), VERSION_PAR)) {
    fprintf(stderr, "\nIncorrect parameter file version : %s\n", getC("EuGene.versionPAR"));
    fprintf(stderr,"Version %s required\n", VERSION_PAR);
    exit(2);
  }
  else
    fprintf(stderr, "(%s)\n", m["EuGene.versionPAR"]);
  
  // Remplir le tableau des longueurs min de chaque etat (+ 6 pour les Single)
  m["EuGene.minL0"] = m["EuGene.minL1"] = m["EuGene.minL2"] = m["EuGene.minEx"];
  m["EuGene.minL3"] = m["EuGene.minL4"] = m["EuGene.minL5"] = m["EuGene.minEx"];
  m["EuGene.minL6"] = m["EuGene.minL7"] = m["EuGene.minL8"] = m["EuGene.minIn"];
  m["EuGene.minL9"] = m["EuGene.minL10"]= m["EuGene.minL11"]= m["EuGene.minIn"];
  m["EuGene.minL12"]= m["EuGene.minL13"]= m["EuGene.minL14"]= m["EuGene.minSg"];
  m["EuGene.minL15"]= m["EuGene.minL16"]= m["EuGene.minL17"]= m["EuGene.minSg"];
}

// ------------------------
//  getChar param.
// ------------------------
char* Parameters :: getC(char *key, int index = 0)
{
  if (!index && m.count(key))
    return (char*)m[key];

  int len = strlen(key);
  char *altkey = new char[len+10];

  strcpy(altkey,key);
  key = altkey+len;
  sprintf(key,"[%d]",index);

  if (m.count(altkey)) {
    key = (char*)m[altkey];
    delete altkey;
    return key;
  }
  else {
    fprintf(stderr,"WARNING: Undefined key %s\n",altkey);
    return NULL;
  }
}

// ------------------------
//  getDouble param.
// ------------------------
double Parameters :: getD(char *key, int index = 0)
{
  char *res = getC(key,index);

  if (res) {
    if(!strcmp(res, "NINFINITY"))
      return NINFINITY;
    else
      return atof(res);
  }
  return 0.0;
}

// ------------------------
//  getInt param.
// ------------------------
int Parameters :: getI(char *key,int index = 0)
{
  char *res = getC(key,index);

  if(res) {
    if(!strcmp(res, "TRUE"))
      return TRUE;
    else if(!strcmp(res, "FALSE"))
      return FALSE;
    else
      return atoi(res);
  }
  return 0; 
}
  
// ------------------------
//  Get Use.Sensor.
// ------------------------
int Parameters :: getUseSensor(char **key, int *val)
{
  int l;
  
  while(iter != m.end()) {
    while(iter != m.end() &&
	  (iter->first[0] != 'S' || iter->first[1] != 'e' ||
	   iter->first[2] != 'n' || iter->first[3] != 's' ||
	   iter->first[4] != 'o' || iter->first[5] != 'r'))
      ++iter;
    if(iter != m.end()) {
      l = strlen(iter->first) - 1;
      if(iter->first[l-2] == 'u' && iter->first[l-1] == 's' && iter->first[l] == 'e')  
	++iter;
      else {
	*key = (char*)iter->first;
	*val = atoi(iter->second);
	++iter;
	return TRUE;
      }
    }
  }
  return FALSE;
}

// ------------------------
//  set param.
// ------------------------
void Parameters :: set (const char *key, const char *value)
{
  m[key] = value;
}

// ------------------------
//  Default destructor.
// ------------------------
Parameters :: ~Parameters ()
{
  m.clear();
}

