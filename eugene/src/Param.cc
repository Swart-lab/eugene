#include "Param.h"

// ------------------------
// Default constructor.
// ------------------------
Parameters :: Parameters ()
{
}

// ------------------------
// InitParam.
// ------------------------
void Parameters :: InitParam (int argc, char * argv[])
{
  ReadArg(argc, argv);

  fprintf(stderr,"EuGene rel. %s\n",VERSION);
  fprintf(stderr,"Loading parameters file...");
  fflush(stderr);

  ReadPar(argv[0]);
}

// ------------------------
// Read arguments line.
// ------------------------
void Parameters :: ReadArg(int argc, char * argv[])
{
  int carg, errflag;
  
  DFT_MATRIX = "default.mat";
  DFT_OUTPUT = "./";
  
  // prior on the initial state, selon Sato et al 1999 / Terryn et
  // al. 1999
  ExonPrior  = 0.33, IntronPrior = 0.17;
  InterPrior = 0.4,  FivePrimePrior = 0.03, ThreePrimePrior = 0.07;
  
  // process args (default values)
  glen = golap = gfrom = gfromSave = gto = gtoSave = -1;
  resx = 900;                 // x res for graph. output
  resy = 400;                 // y res for graph. output
  graph    = FALSE;           // don't produce a graphical output
  estopt   = FALSE;           // don't try to read a EST file
  estanal  = FALSE;           // don't try to analyze EST support
  userinfo = FALSE;           // shall we read a user info file
  raflopt  = FALSE;           // don't try to read a est.rafl file
  blastopt = FALSE;           // don't try to read a blast file
  ncopt    = FALSE;           // don't try to read a non coding input
  normopt = 1;                // normalize across frames
  window  = 97;               // window length
  printopt = 'l';             // short print format 
  offset  = 0;                // no offset
  
  (void) strcpy(matname, DFT_MATRIX); // default matrix file
  (void) strcpy(outputname, DFT_OUTPUT); // default output
  *fstname = '\000';                  // no default input
  parname[0] = '0';                   // no -P
  errflag = 0;
  
  while ((carg = getopt(argc, argv, "UREdrshm:w:f:n:o:p:x:y:c:u:v:g::b::l:P:O:")) != -1) {
    
    switch (carg) {
      
    case 'n':           /* -n normalize across frames      */
      if (! GetIArg(optarg, &normopt, normopt))
	errflag++;
      break;
      
    case 'h':           /* help                             */
      errflag++;
      break;
      
    case 's':           /* Single gene mode. Start/End on IG only   */
      ExonPrior = 0.0;
      IntronPrior = 0.0;
      //      InterPrior = 1.0;
      FivePrimePrior = 0.0;
      ThreePrimePrior = 0.0;
      break;

    case 'r':    /* RepeatMasker input */
      ncopt = TRUE;
      break;
      
    case 'c':           /* -c couverture      */
      if (! GetIArg(optarg, &golap, golap))
        errflag++;
      break;
      
    case 'l':           /* -l imglen      */
      if (! GetIArg(optarg, &glen, glen))
        errflag++;
      break;
      
    case 'u':           /* -u From      */
      if (! GetIArg(optarg, &gfrom, gfrom))
	errflag++;
      gfromSave = gfrom;
      break;
      
    case 'v':           /* -v To      */
      if (! GetIArg(optarg, &gto, gto))
	errflag++;
      gtoSave = gto;
      break;
      
    case 'x':           /* -x resx      */
      if (! GetIArg(optarg, &resx, resx))
        errflag++;
      break;
      
    case 'y':           /* -y resy      */
      if (! GetIArg(optarg, &resy, resy))
        errflag++;
      break;
      
    case 'g':          /* -g "Graphic File"                */
      if (optarg) {
	if (strspn(".fasta", optarg)    != 6 
	    && strspn(".fsa", optarg)   != 4
	    && strspn(".tfa", optarg)   != 4
	    && strspn(".txt", optarg)   != 4
	    && optarg[0] != '-') {
	  if(strcmp(optarg, "NO"))      // Argument suivant le -g
	    strcpy(grnameArg,optarg);   // != NO pris en compte
	  else 
	    grnameArg[0] = 0;           // == NO non pris en compte
	  graph = TRUE;
	}
	else errflag++;
      }
      else
	grnameArg[0] = 0;
      break;
      
    case 'p':          /* print opt: short/long/detailed   */
      printopt = optarg[0];
      if ((printopt == 'h') &&  graph == FALSE) {
	graph = TRUE; // HTML output means graphical output 
	grnameArg[0] = 0;
      }
      if ((printopt != 's') && (printopt != 'l') && (printopt != 'g') && 
	  (printopt != 'd') && (printopt != 'h'))
	errflag++;
      break;
      
    case 'm':           /* -m matrix                        */
      (void) strcpy(matname, optarg);
      break;

    case 'O':           /* -O output                        */
      (void) strcpy(outputname, optarg);
      break;

    case 'P':           /* -P .par                          */
      (void) strcpy(parname, optarg);
      break;

    case 'o':           /* -o offset                        */
      if (! GetIArg(optarg, &offset, offset))
	errflag++;
      break;
      
    case 'w':           /* -w window                        */
      if (! GetIArg(optarg, &window, window))
	errflag++;
      if (window == 0) window = 1;
      break;

    case 'b':           /* -b  use blastx result            */
      if (optarg) {
	(void) strcpy(blastArg, optarg);
	blastopt = TRUE;
	if(strlen(blastArg) > 8)
	  errflag++;
      }
      else {
	blastArg[0] = '1';
	blastArg[1] = '2';
	blastArg[2] = '3';
      }
      break;
      
    case 'E':
      estanal = TRUE;
      break;
      
    case 'U':
      userinfo = TRUE;
      break;

    case 'd':           /* -d  use cDNA blastn results      */
      estopt = TRUE;
      break; 

    case 'R':           /* -R use RAFL-like EST*/
      raflopt = TRUE;
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
// Read parameters file.
// ------------------------
void Parameters :: ReadPar(char *argv)
{
  char Line    [MAX_LINE];
  char tempname[FILENAME_MAX+1];
  int  nbPar = 0;

  if(parname[0] == '0') {
    strcpy(tempname, argv);
    strcat(tempname, ".par");
    fp = FileOpen(EugDir, BaseName(tempname),"r");
  }
  else {
    fp = FileOpen(EugDir, parname,"r");
    strcpy(tempname, parname);
  }

  fgets (Line, MAX_LINE, fp);
  fgets (Line, MAX_LINE, fp);
  fgets (Line, MAX_LINE, fp);
  if (fgets (Line, MAX_LINE, fp) != NULL)
  if (sscanf(Line, "%s", clef) == 1)
    if (fgets (Line, MAX_LINE, fp) != NULL)
    if (sscanf(Line, "%s", versionPAR) == 1)
      if (fgets (Line, MAX_LINE, fp) != NULL)
      if (sscanf(Line, "%lf", &FsP) == 1)
        if (fgets (Line, MAX_LINE, fp) != NULL)
        if (sscanf(Line, "%lf %lf", &StartP, &StartB) == 2)
          if (fgets (Line, MAX_LINE, fp) != NULL)
          if (sscanf(Line, "%lf", &StopP) == 1)
            if (fgets (Line, MAX_LINE, fp) != NULL)
            if (sscanf(Line, "%lf", &TransStartP) == 1)
              if (fgets (Line, MAX_LINE, fp) != NULL)
              if (sscanf(Line, "%lf", &TransStopP) == 1)
       	        if (fgets (Line, MAX_LINE, fp) != NULL)
                if (sscanf(Line, "%lf %lf", &AccP[0], &AccB[0]) == 2)
                  if (fgets (Line, MAX_LINE, fp) != NULL)
     	          if (sscanf(Line, "%lf %lf", &DonP[0], &DonB[0]) == 2)
                    if (fgets (Line, MAX_LINE, fp) != NULL)
    	            if (sscanf(Line, "%lf %lf", &AccP[1], &AccB[1]) == 2)
                      if (fgets (Line, MAX_LINE, fp) != NULL)
    	              if (sscanf(Line, "%lf %lf", &DonP[1], &DonB[1]) == 2)
			if (fgets (Line, MAX_LINE, fp) != NULL)
			if (sscanf(Line, "%lf", &EstP) == 1)
			  if (fgets (Line, MAX_LINE, fp) != NULL)
			  if (sscanf(Line, "%d", &EstM) == 1)
			    if (fgets (Line, MAX_LINE, fp) != NULL)
			    if (fgets (Line, MAX_LINE, fp) != NULL)
			    if (sscanf(Line, "%lf %lf %lf %lf %lf %lf %lf %lf",
				       &BlastS[0],&BlastS[1],&BlastS[2],&BlastS[3],
				       &BlastS[4],&BlastS[5],&BlastS[6],&BlastS[7]) == 8)
			      if (fgets (Line, MAX_LINE, fp) != NULL)
			      if (fgets (Line, MAX_LINE, fp) != NULL)
			      if (sscanf(Line, "%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d",
					 &MinEx, &MinIn, &MinSg, &MinFlow, &MinConv,
					 &MinDiv, &MinFivePrime, &MinThreePrime) == 8)
				if (fgets (Line, MAX_LINE, fp) != NULL)
				if (fgets (Line, MAX_LINE, fp) != NULL)
				if (sscanf(Line, "%lf %lf %lf %lf %lf",
					   &TransitionWeight[0], &TransitionWeight[1],
					   &TransitionWeight[2], &TransitionWeight[3],
					   &TransitionWeight[4]) == 5) {
				  nbPar = TRUE;
				  fprintf(stderr,"done\n");
				  if(strcmp(versionPAR,VERSION_PAR)) {
				    fprintf(stderr,"Incorrect parameter file version : %s\n",
					    versionPAR);
				    fprintf(stderr,"Version %s required\n", VERSION_PAR);
				    exit(2);
				  }
				  else
				    fprintf(stderr,"Parameters file %s\n", versionPAR);
				}
  if (!nbPar) {
    if (EugDir != NULL)
      fprintf(stderr, "Incorrect parameter file %s/%s\n",EugDir,BaseName(tempname));
    else
      fprintf(stderr, "Incorrect parameter file %s\n",BaseName(tempname));
    exit(2);
  }

  // Remplir le tableau des longueurs min de chaque etat (+ 6 pour les Single)
  for (int i = 0; i < 6; i++) {
    MinLength[i]    = MinEx;
    MinLength[i+6]  = MinIn;
    MinLength[i+12] = MinSg;
  }
}

// ------------------------
//  Default destructor.
// ------------------------
Parameters :: ~Parameters ()
{
}
