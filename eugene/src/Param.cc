#include "Param.h"

// ------------------------
//  Default constructor.
// ------------------------
Parameters :: Parameters ()
{
  DFT_MATRIX = "default.mat";
  DFT_OUTPUT = "./";
}

// ------------------------
//  initParam.
// ------------------------
void Parameters :: initParam (int argc, char * argv[])
{
  ReadArg(argc, argv);

  fprintf(stderr,"EuGene rel. %s\n",VERSION);
  fprintf(stderr,"Loading parameters file...");
  fflush(stderr);

  ReadPar(argv[0]);
 
  //for (iter = m.begin(); iter!=m.end(); ++iter)
  //fprintf(stderr,"%s  est associé à %s\n",iter->first, iter->second);

  iter = m.begin();  // Cf. : getUseSensor
}

// ------------------------
//  Read arguments line.
// ------------------------
void Parameters :: ReadArg(int argc, char * argv[])
{
  int carg, errflag;
  
  // Prior on the initial state, Sato et al 1999 / Terryn et al. 1999
  m["ExonPrior"]   = "0.33", m["IntronPrior"]    = "0.17";
  m["InterPrior"]  = "0.4",  m["FivePrimePrior"] = "0.03";
  m["ThreePrimePrior"] = "0.07";
  
  // Process args (default values)
  m["glen"]      = m["golap"] = m["gfrom"]   = "-1";
  m["gfromSave"] = m["gto"]   = m["gtoSave"] = "-1";
  m["resx"]     = "900";        // x res for graph. output
  m["resy"]     = "400";        // y res for graph. output
  m["graph"]    = "FALSE";      // don't produce a graphical output
  m["estopt"]   = "FALSE";      // don't try to read a EST file
  m["estanal"]  = "FALSE";      // don't try to analyze EST support
  m["userinfo"] = "FALSE";      // shall we read a user info file
  m["raflopt"]  = "FALSE";      // don't try to read a est.rafl file
  m["blastopt"] = "FALSE";      // don't try to read a blast file
  m["ncopt"]    = "FALSE";      // don't try to read a non coding input
  m["normopt"]  = "1";          // normalize across frames
  m["window"]   = "97";         // window length
  m["printopt"] = "l";          // short print format 
  m["offset"]   = "0";          // no offset
  m["matname"]    = DFT_MATRIX; // default matrix file
  m["outputname"] = DFT_OUTPUT; // default output
  m["fstname"] = "\000";        // no default input
  parname[0]   = '0';           // no -P
  FsP          = 0;
  errflag      = 0;

  while ((carg = getopt(argc, argv, "UREdrshm:w:f:n:o:p:x:y:c:u:v:g::b::l:P:O:")) != -1) {
    
    switch (carg) {
      
    case 'n':           /* -n normalize across frames      */
      if (! GetIArg(optarg, &normopt, normopt))
	errflag++;
      else m["normopt"] = intToChar(normopt);
      break;
      
    case 'h':           /* help                             */
      errflag++;
      break;
      
    case 's':           /* Single gene mode. Start/End on IG only   */
      m["ExonPrior"]   = "0.0";
      m["IntronPrior"] = "0.0";
      //InterPrior = 1.0;
      m["FivePrimePrior"]  = "0.0";
      m["ThreePrimePrior"] = "0.0";
      break;

    case 'r':    /* RepeatMasker input */
      m["ncopt"] = "TRUE";
      break;
      
    case 'c':           /* -c couverture      */
      if (! GetIArg(optarg, &golap, golap))
        errflag++;
      else m["golap"] = intToChar(golap);
      break;
      
    case 'l':           /* -l imglen      */
      if (! GetIArg(optarg, &glen, glen))
        errflag++;
      else m["glen"] = intToChar(glen);
      break;
      
    case 'u':           /* -u From      */
      if (! GetIArg(optarg, &gfrom, gfrom))
	errflag++;
      else m["gfrom"] = intToChar(gfrom);
      m["gfromSave"] = m["gfrom"];
      break;
      
    case 'v':           /* -v To      */
      if (! GetIArg(optarg, &gto, gto))
	errflag++;
      else m["gto"] = intToChar(gto);
      m["gtoSave"] = m["gto"];
      break;
      
    case 'x':           /* -x resx      */
      if (! GetIArg(optarg, &resx, resx))
        errflag++;
      else m["resx"] = intToChar(resx);
      break;
      
    case 'y':           /* -y resy      */
      if (! GetIArg(optarg, &resy, resy))
        errflag++;
      else m["resy"] = intToChar(resy);
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
	  m["graph"] = "TRUE";
	}
	else errflag++;
      }
      else
	m["grnameArg"] = "0";
      break;
      
    case 'p':           /* print opt: short/long/detailed   */
      m["printopt"] = optarg;
      if ((m["printopt"] == "h") &&  m["graph"] == "FALSE") {
	m["graph"]     = "TRUE"; // HTML output means graphical output 
	m["grnameArg"] = "0";
      }
      if ((m["printopt"][0] != 's') && (m["printopt"][0] != 'l') &&
	  (m["printopt"][0] != 'g') && (m["printopt"][0] != 'd') &&
	  (m["printopt"][0] != 'h'))
	errflag++;
      break;
      
    case 'm':           /* -m matrix                        */
      m["matname"] = optarg;
      break;

    case 'O':           /* -O output                        */
      m["outputname"] = optarg;
      break;

    case 'P':           /* -P .par                          */
      (void) strcpy(parname, optarg);
      break;

    case 'o':           /* -o offset                        */
      if (! GetIArg(optarg, &offset, offset))
	errflag++;
      else m["offset"] = intToChar(offset);
      break;
      
    case 'w':           /* -w window                        */
      if (! GetIArg(optarg, &window, window))
	errflag++;
      else m["window"] = intToChar(window);
      if (m["window"] == "0") m["window"] = "1";
      break;

    case 'b':           /* -b  use blastx result            */
      if (optarg) {
	m["blastArg"] = optarg;
	m["blastopt"] = "TRUE";
	if(strlen(m["blastArg"]) > 8)
	  errflag++;
      }
      else m["blastArg"] = "123";
      break;
      
    case 'E':
      m["estanal"] = "TRUE";
      break;
      
    case 'U':
      m["userinfo"] = "TRUE";
      break;

    case 'd':           /* -d  use cDNA blastn results      */
      m["estopt"] = "TRUE";
      break; 

    case 'R':           /* -R use RAFL-like EST*/
      m["raflopt"] = "TRUE";
      break;

    case 'f':           /* -f frameshift loglike            */
      if (! GetDArg(optarg, &FsP, FsP))
	errflag++;
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
    fprintf(stderr, "Usage: EuGene [-h] [-m matrix] [-P .par] [-n 0|1|2] [-s] [-p h|g|s|l|d]\n"
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
  
  if(parname[0] == '0') {
    strcpy(tempname, argv);
    strcat(tempname, ".par");
    fp = FileOpen(EugDir, BaseName(tempname),"r");
  }
  else {
    fp = FileOpen(EugDir, parname,"r");
    strcpy(tempname, parname);
  }

  fgets (line, MAX_LINE, fp);
  fgets (line, MAX_LINE, fp);
  fgets (line, MAX_LINE, fp);

  while(fgets (line, MAX_LINE, fp) != NULL) {
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

  fclose(fp);
  fprintf(stderr,"done\n");
  if(strcmp(getC("EuGene.versionPAR"), VERSION_PAR)) {
    fprintf(stderr, "Incorrect parameter file version : %s.\n", getC("EuGene.versionPAR"));
    fprintf(stderr,"Version %s required\n", VERSION_PAR);
    exit(2);
  }
  else
    fprintf(stderr, "Parameters file %s\n", m["EuGene.versionPAR"]);
  
  // Remplir le tableau des longueurs min de chaque etat (+ 6 pour les Single)
  m["EuGene.minL0"] = m["EuGene.minL1"] = m["EuGene.minL2"] = m["EuGene.minEx"];
  m["EuGene.minL3"] = m["EuGene.minL4"] = m["EuGene.minL5"] = m["EuGene.minEx"];
  m["EuGene.minL6"] = m["EuGene.minL7"] = m["EuGene.minL8"] = m["EuGene.minIn"];
  m["EuGene.minL9"] = m["EuGene.minL10"]= m["EuGene.minL11"]= m["EuGene.minIn"];
  m["EuGene.minL12"]= m["EuGene.minL13"]= m["EuGene.minL14"]= m["EuGene.minSg"];
  m["EuGene.minL15"]= m["EuGene.minL16"]= m["EuGene.minL17"]= m["EuGene.minSg"];

  // La ligne d'argument est prioritaire sur le .par
  // Si info défini en argument on modifie sa valeur lu dans .par
  if( FsP != 0 )                                 // -f -> Frameshift
    m["EuGene.frameshift"] = doubleToChar(FsP);
  if( m["estopt"] == "TRUE" )                    // -d -> Est
    m["Sensor.Est.use"] = "TRUE";
  if( m["blastopt"] == "TRUE" )                  // -b -> BlastX
    m["Sensor.BlastX.use"] = "TRUE";
  if( m["userinfo"] == "TRUE" )                  // -U -> User
    m["Sensor.User.use"] = "TRUE";
  if( m["raflopt"] == "TRUE" )                   // -R -> Riken
    m["Sensor.Riken.use"] = "TRUE";
  if( m["ncopt"] == "TRUE" )                     // -r -> Repeat
    m["Sensor.Repeat.use"] = "TRUE";
}

// ------------------------
//  getDouble param.
// ------------------------
char* Parameters :: getC(char *key)
{
  if(m.count(key))
    return (char*)m[key];
  else {
    fprintf(stderr,"WARNING: Undefined key %s\n",key);
    return (char*)"";
  }
}

// ------------------------
//  getDouble param.
// ------------------------
double Parameters :: getD(char *key)
{
  if(m.count(key)) {
    if(!strcmp(m[key], "NINFINITY"))
      return NINFINITY;
    else
      return atof(m[key]);
  }
  else {
    fprintf(stderr,"WARNING: Undefined key %s\n",key);
    return 0.0;
  }
}

// ------------------------
//  getInt param.
// ------------------------
int Parameters :: getI(char *key)
{
  if(m.count(key)) {
    if(!strcmp(m[key], "TRUE"))
      return TRUE;
    else if(!strcmp(m[key], "FALSE"))
      return FALSE;
    else
      return atoi(m[key]);
  }

  else {
    fprintf(stderr,"WARNING: Undefined key %s\n",key);
    return 0; 
  }
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
//  Convert int -> char*.
// ------------------------
const char* Parameters :: intToChar(int val)
{
  strstream buf;
  string s;
  buf << val;
  buf >> s;  
  return s.c_str();
}

// --------------------------
//  Convert double -> char*.
// --------------------------
const char* Parameters :: doubleToChar(double val)
{
  strstream buf;
  string s;
  buf << val;
  buf >> s;  
  return s.c_str();
}

// ------------------------
//  Default destructor.
// ------------------------
Parameters :: ~Parameters ()
{
  m.clear();
}
