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
  fprintf(stderr,"EuGene rel. %s",VERSION);
  fflush(stderr);

  ReadPar(argv[0]);
  ReadArg(argc, argv);

  if (m.count("Param.debug")) 
    for (iter = m.begin(); iter!=m.end(); ++iter)
      fprintf(stderr,"%s = %s\n",iter->first, iter->second);

  iter = m.begin();  // Cf. : getUseSensor
  
  fprintf(stderr,"-%s (%s)\n",getC("EuGene.organism"),VERSION_DATE);
  fprintf(stderr,"Parameters file loaded.\n");
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

  while ((carg = getopt(argc, argv, "GUREBdrshm:w:f:n:o:p:x:y:c:u:v:g::b::l:O:D:t::M:Z")) != -1) {
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
	fprintf(stderr,"\nInvalid command line parameter definition: %s\n",optarg);
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

    case 'O':           /* -O output                        */
      m["Output.Prefix"] = optarg;
      break;
      
    case 'p':           /* -p print opt: s -> short, l -> long, d -> detailed
                                         h -> html,  H -> HTML(EuGeneHom)
                                         g -> GFF,   a -> araset  */
      m["Output.format"] = optarg;
      if ((optarg[0] == 'h') && getI("Output.graph") == 0) {
	m["Output.graph"] = "TRUE"; // HTML output means graphical output 
	m["grnameArg"]    = "0";
      }
      if ((optarg[0] != 's') && (optarg[0] != 'l') && (optarg[0] != 'g') &&
	  (optarg[0] != 'd') && (optarg[0] != 'h') && (optarg[0] != 'a') &&
          (optarg[0] != 'H'))
	errflag++;
      break;
      
    case 'm':           /* -m IMM matrix                    */
      m["Sensor.MarkovIMM.use"] = "TRUE" ;
      m["MarkovIMM.matname"] = optarg;
      break;

    case 'M':           /* -M proteic markovian matrix      */
      m["Sensor.MarkovProt.use"] = "TRUE";
      m["MarkovProt.matname"] = optarg;
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
	if(strlen(optarg) > 10)
	  errflag++;
      }
      break;

    case 'E':
      m["Est.PostProcess"] = "TRUE";
      break;

    case 'B':
      m["BlastX.PostProcess"] = "TRUE";
      break;

    case 'U':
      m["Sensor.User.use"] = "TRUE";
      break;

    case 'G':           /* -G use sensor GFF */
      m["Sensor.GFF.use"] = "TRUE";
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

    case 'Z':           // parameters optimization
      m["ParaOptimization.Use"] = "TRUE";
      if (optarg) m["ParaOptimization.TrueCoordFile"] = optarg;
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
    fprintf(stderr, "\nUsage: EuGene [-h] [-m matrix] [-n 0|1|2] [-s] [-p h|g|s|l|d|a]\n"
	    "              [-w window] [-b {levels}] [-d] [-R] [-E] [-U] [-o offset]\n"
	    "              [-g {graphArg}] [-G] [-u start] [-v end] [-l len] [-c olap]\n"
	    "              [-x xres] [-y yres] [-Z coord] FASTA files\n");
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
  int  n=0;

  strcpy(tempname, argv);
  strcat(tempname, ".par");
  fp = FileOpen(EugDir, BaseName(tempname),"r");

  while(fgets (line, MAX_LINE, fp) != NULL) {
    n++;
    if (line[0] != '#') {
      key = new char[FILENAME_MAX+1];
      val = new char[FILENAME_MAX+1];
      if(sscanf(line, "%s %s", key, val) == 2)  
	m[key] = val;
      else {
	if (EugDir != NULL)
	  fprintf(stderr, "\nIncorrect parameter file %s/%s\n",EugDir,BaseName(tempname));
	else
	  fprintf(stderr, "\nIncorrect parameter file %s line %d\n",BaseName(tempname),n);
	exit(2);
      }
    }
  }

  fclose(fp);

  if(strcmp(getC("EuGene.versionPAR"), VERSION_PAR)) {
    fprintf(stderr, "\nIncorrect parameter file version : %s\n", getC("EuGene.versionPAR"));
    fprintf(stderr,"Version %s required\n", VERSION_PAR);
    exit(2);
  }
  
  // Remplir le tableau des longueurs min de chaque etat (+ 6 pour les Single)
  m["EuGene.minL0"] = m["EuGene.minL1"] = m["EuGene.minL2"] = m["EuGene.minEx"];
  m["EuGene.minL3"] = m["EuGene.minL4"] = m["EuGene.minL5"] = m["EuGene.minEx"];
  m["EuGene.minL6"] = m["EuGene.minL7"] = m["EuGene.minL8"] = m["EuGene.minIn"];
  m["EuGene.minL9"] = m["EuGene.minL10"]= m["EuGene.minL11"]= m["EuGene.minIn"];
  m["EuGene.minL12"]= m["EuGene.minL13"]= m["EuGene.minL14"]= m["EuGene.minSg"];
  m["EuGene.minL15"]= m["EuGene.minL16"]= m["EuGene.minL17"]= m["EuGene.minSg"];
}

// ------------------------
//  count.
// ------------------------
int Parameters :: count(char *key)
{
  return ( m.count(key) );
}

// ------------------------
//  getChar param.
// ------------------------
char* Parameters :: getC(char *key, int index)
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
    delete [] altkey;
    return key;
  }
  else {
    fprintf(stderr,"\nError: Undefined key %s\n",altkey);
    exit(2);
  }
}

// ------------------------
//  getDouble param.
// ------------------------
double Parameters :: getD(char *key, int index)
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
int Parameters :: getI(char *key,int index)
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
  char *s; s = new char[FILENAME_MAX+1];

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
	strcpy(s,iter->first); strcat(s,".use"); 
	if (!strcmp(m[s],"TRUE")) {
	  *key = (char*)iter->first;
	  *val = atoi(iter->second);
	  ++iter;
	  return TRUE;
	}
	++iter;
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

void Parameters :: setD (const char *key, double n)
{
  char *buffer = new char[FILENAME_MAX+1];
  sprintf (buffer, "%10f",n);
  m[key] = buffer;
}




// ------------------------
//  Default destructor.
// ------------------------
Parameters :: ~Parameters ()
{
  m.clear();
}

// ----------------------------------------------------------
// WriteParam: write a new parameter file, named <para_file>.<date>.OPTI.par,
//             modifying the value of parameters in para_name
//             with the new value in para_val
// BEWARE: the comments after a parameter value change are not kept
//         in case of no comment the result of strchr(line,'#') was: (null).
// Evaluation: the name of the written parameter file 
// -----------------------------------------------------------
std::string Parameters::WriteParam (const char* para_file, std::vector<std::string> para_name, 
				    std::vector<double> para_val)
{
  FILE   *fp, *fp_opti;
  char line[MAX_LINE], new_line[MAX_LINE];
  char key[MAX_LINE], val[MAX_LINE];
  char filename[FILENAME_MAX+1];
  unsigned int i;
  bool find_para;
  char *d = new char[MAX_LINE];

  strcpy(filename, para_file);
  strcat(filename, ".par");
  fp = FileOpen(EugDir, BaseName(filename),"r");

  strcpy(filename, para_file);
  strcat(filename, ".");
  GetStrDate(d);
  strcat(filename, d);
  strcat(filename, ".OPTI");
  strcat(filename, ".par");
  fp_opti = FileOpen(EugDir, BaseName(filename),"w");
  
  while(fgets (line, MAX_LINE, fp) != NULL) {
    find_para = false;
    if (line[0] != '#') 
      if (sscanf(line, "%s %s", key, val)==2) 
	for (i=0; i<para_name.size(); i++) 
	  if ( para_name[i] == (std::string) key ) {
	    find_para = true;
	    sprintf(new_line,"%s\t%f\n", key, para_val[i]);
	    i = para_name.size();
	  }
    if (!find_para)
      for (int i=0; i<MAX_LINE; i++) new_line[i] = line[i];
    fprintf(fp_opti,"%s",new_line);
  }

  fclose(fp);
  fclose(fp_opti);

  return (std::string) filename;
}



// -----------------------------------------------------------
// -----------------------------------------------------------
void Parameters::ResetIter(void)
{
  iter = m.begin();
}
