//  This program finds exons/introns and intergenic regions (including UTR)
//  Copyright T. Schiex 1999

// Wed Mar  3 12:06:19 1999 : des starts factices en extremites de seq
// Wed Mar  3 12:28:08 1999 : Indep_eval decale de 1 
// Mon Mar  8 13:27:04 1999 : abaisse la force des START (scale)
// Mon Mar  8 17:27:12 1999 : supprime s/ es. Data commence a 0
// Mon Mar 15 10:12:30 1999 : passage a une version eucaryote
// Mon Mar 15 15:37:00 1999 : patch sigsegv annick - delete Data => (Data-1)
// Thu Mar 18 23:01:35 1999 : usage de macs decale
// Thu Mar 18 23:03:35 1999 : prise en compte de HITS BLAST EST/PROTEINES
// Sun May 23 15:30:26 1999 : modele intronique
// Tue Jun 22 16:50:21 1999 : modele intronique trop fort. 
//                            Utilise en intergenique too
// Tue Jul 6  13:01:23 1999 : modele intronique dirige,
//                            modele intergenique non.
// Mon Jul 19 16:07:54 1999 : introduction de StartB et separation
//                            des parametres Acc et Don
// Mon Jul 19 16:09:54 1999 : Le Start est traite comme une proba cond.
// Mon Oct 18 14:58:36 1999 : macs => backP, traitement des longueurs min


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <ctype.h>
#include <limits.h>
#include <float.h>
#include <time.h>
#include <assert.h>
#include <errno.h>
#include <getopt.h>

//#define STAND 1

#define DFT_MATRIX          "default.mat"
#define REAL double

#ifndef FILENAME_MAX
#define FILENAME_MAX        1024
#endif

#define NORM(x,n) (((n)+(MAX(-(n),x)))/(n))


//#define SpliceNOnly

// ------------------ Globals ---------------------
const int  MODEL_LEN = 9;
const int  SIMPLE_MODEL_LEN = 6;
const int  ALPHABET_SIZE = 4;
const REAL NINFINITY = log(0.0);
const int PredWidth = 2;

// Les Hits EST

const char HitForward    = 0x1;
const char MLeftForward  = 0x2;
const char GapForward    = 0x4;
const char MRightForward = 0x8; 

const char HitReverse    = 0x10;
const char MLeftReverse  = 0x20;
const char GapReverse    = 0x40;
const char MRightReverse = 0x80; 

const char Hit           = HitForward    | HitReverse;
const char MLeft         = MLeftForward  | MLeftReverse;
const char Gap           = GapForward    | GapReverse;
const char MRight        = MRightForward | MRightReverse;

const char NotAHit       = MLeft | MRight | Gap;

#define Inconsistent(x) (((x) & Hit) && ((x) & Gap))

const short int State2Phase[13] = {1,2,3,-1,-2,-3,4,4,4,-4,-4,-4,0};

REAL Score [9],NScore[9];
long int Data_Len;
int normopt,blastopt,estopt;
int window, offset,graph,resx,resy;
int    gfrom,gto,golap,glen;

extern char     *optarg;   
extern int      optind;

#include "System.h"
#include "gene.h"
#include "Stop.h"
#include "BStrArray.h"
#include "EuIMMScoreAS.h"
#include "ArgDecode.h"
#include "EuStart.h"
#include "BackP.h"
#ifdef SplicePOnly
#include "SpliceP.h"
#else
#ifdef SpliceNOnly
#include "SpliceN.h"
#else 
#include "SpliceNP.h"
#endif
#endif
#include "../GDIF/gdIF.h"

// Les longueurs. Au pif (a part les introns)

const int MinIgFlow = 250;
const int MinIgConv = 150;
const int MinIgDiv  = 300;

const int MinEx = 0;
const int MinIn = 50;
const int MinSg = 150;

// tableau des longueurs min de chaque etat (+ 6 pour les Single)
const int MinLength[19] = 
   {MinEx, MinEx, MinEx,MinEx, MinEx, MinEx,
    MinIn, MinIn, MinIn,MinIn, MinIn, MinIn,
    MinIgFlow, 
    MinSg, MinSg,MinSg, MinSg, MinSg,MinSg};

const unsigned char SwitchMask[19] = 
   {SwitchAny,SwitchAny,SwitchAny,SwitchAny,SwitchAny,SwitchAny,
    SwitchAny,SwitchAny,SwitchAny,SwitchAny,SwitchAny,SwitchAny,
    SwitchAny,
    SwitchStart,SwitchStart,SwitchStart,SwitchStop,SwitchStop,SwitchStop,
    };

// les meme en strict (pour les intergeniques)
const int StrictMinLength[4] =
   {MinIgFlow, MinIgDiv, MinIgFlow, MinIgConv};

const unsigned char StrictSwitchMask[4] = 
   {SwitchStop,SwitchStart,SwitchStart,SwitchStop};


// ------------------ Globals ---------------------

int IsPhaseOn(char p, int pp)
{
  if (p == 6) return FALSE;
  if (p == pp) return TRUE;
  return FALSE;
}


// Convertit les phases 0-6 en 1 2 3 -1 -2 -3 0
int PhaseAdapt(char p)
{
  if (p == 12) return 0;
  else if (p < 3) return (1+p);
  else if (p < 6) return (2-p);
  else if (p < 9) return (p-2);
  else return (5-p);
  
}

void PrintPhase(char p)
{
  printf("%2d ", PhaseAdapt(p));
  return;
}

// Convertit les phases 1 2 3 -1 -2 -3 0 en 0-6
char ph06(char p)
{
  if (p == 0) return 6;
  else if (p > 0) return (p-1);
  else return 2-p;   
}
// -------------------------------------------------------------------------
// sortie graphique
// -------------------------------------------------------------------------
void PlotPredictions(char *Choice,REAL *BaseScore[9], 
		     unsigned char *Stop[2],REAL *Start[2],
		     REAL *Acc[2], REAL *Don[2])
{
  static double LScore[9] = {-1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0,-1.0,-1.0};  
  int i,j,p;

  for (i=0; i<Data_Len; i++)
    PlotBarI(i,State2Phase[Choice[i]],0.4,PredWidth,1);
  
  for (j = 0 ; j < 9 ; j++)
    NScore [j] = 0.0;
  
  for (j = 0 ; j < 9 ; j++)
    for (i = 0; i < window/2; i++)
      NScore [j] += log(BaseScore[j][i]);

  for (i=0; i<Data_Len; i++)
    {
      if (i-window/2 >= 0)
	for (j = 0 ; j < 9 ; j++)
	  NScore [j] -= log(BaseScore[j][i-window/2]);
      
      if (i+window/2 < Data_Len)
	for (j = 0 ; j < 9 ; j++)
	  NScore [j] += log(BaseScore[j][i+window/2]);
      
      for (j = 0 ; j < 9 ; j++) Score[j] = NScore[j];
      
      AmplifyScore(Score,normopt);
      
      if (LScore[0] < 0) for (j = 0 ; j < 9 ; j++) LScore[j] = Score[j];

      p = ((i == 0) ? 0 : i-1);
      
      PlotLine(p,i, 1, 1,LScore[0],Score[0],3);
      PlotLine(p,i, 2, 2,LScore[1],Score[1],3);
      PlotLine(p,i, 3, 3,LScore[2],Score[2],3);
      PlotLine(p,i,-1,-1,LScore[3],Score[3],3);
      PlotLine(p,i,-2,-2,LScore[4],Score[4],3);
      PlotLine(p,i,-3,-3,LScore[5],Score[5],3);
      PlotLine(p,i, 4, 4,LScore[6],Score[6],3);
      PlotLine(p,i,-4,-4,LScore[7],Score[7],3);
      PlotLine(p,i, 0, 0,LScore[8],Score[8],3);
      
      for (j = 0 ; j < 9 ; j++) LScore[j] = Score[j];
    }

  for (i=0; i<Data_Len; i++)
    {
      if (Stop[0][i]) 
	PlotBarF(i,(i%3)+1,0.1,0.2,1);
      
      if (Stop[1][i+1]) 
	PlotBarF(i,-((Data_Len-i-1)%3)-1,0.1,0.2,1);
      
      if (Start[0][i] > 0.0) 
	PlotBarF(i,(i%3)+1,0.5,NORM(log(Start[0][i]),4),2);
      
      if (Start[1][i+1] > 0.0) 
	PlotBarF(i,-((Data_Len-i-1)%3)-1,0.5,NORM(log(Start[1][i+1]),4),2);
      
      if (Acc[0][i] > 0.0) 
	PlotBarF(i,4,0.5,NORM(log(Acc[0][i]),20),4);
      
      if (Acc[1][i+1] > 0.0) 
	PlotBarF(i,-4,0.5,NORM(log(Acc[1][i+1]),20),4);

      if (Don[0][i] > 0.0) 
	PlotBarF(i,4,0.5,NORM(log(Don[0][i]),20),5);
      
      if (Don[1][i+1] > 0.0) 
	PlotBarF(i,-4,0.5,NORM(log(Don[1][i+1]),20),5);
    }    
}

// -------------------------------------------------------------------------
//            MAIN
// -------------------------------------------------------------------------
int main  (int argc, char * argv [])
{
  int              i, j, k, carg, errflag;
  FILE             *fp;
  REAL             *BaseScore[9],*Start[2] ;
  REAL             *Don[2],*Acc[2];
  unsigned char    *Stop[2];
  unsigned char    *ESTMatch;
  long int         Input_Size;
  REAL             Ch_Ct[ALPHABET_SIZE] = {0.0};
  BString_Array    *IMMatrix[5];
  char             *Data;
  char             printopt;
  char             matname[FILENAME_MAX+1], tempname[FILENAME_MAX+1], fstname[FILENAME_MAX+1];
  char             grname[FILENAME_MAX+1];
  char             seqname[MAX_LINE];
  int              EstM;
  REAL             FsP,StartP,StartB,StopP;
  REAL             AccP,AccB,DonP,DonB,BlastS,EstP;
  REAL             ExonPrior,IntronPrior,InterPrior;
  // process args 
  // default values

  glen = -1;
  golap = -1;
  gfrom = -1;
  gto = -1;
  resx = 900;                 // x res for graph. output
  resy = 400;                 // y res for graph. output
  graph = FALSE;              // don't produce a graphical output
  estopt = FALSE;             // don't try to read a EST file
  blastopt = FALSE;           // don't try to read a blast file
  normopt = 1;                // normalize across frames
  window  = 97;               // window length
  printopt = 'l';             // short print format 
  offset = 0;                 // no offset
  ExonPrior = 0.33;
  IntronPrior = 0.17;
  InterPrior = 0.5;

  (void) strcpy(matname, DFT_MATRIX); // default matrix file
  *fstname = '\000';                  // no default input    
  errflag = 0;

  if (! (fp = OpenFile("EuGene.par","r")))
    {
      fprintf(stderr, "cannot open parameter file EuGene.par\n");
      exit(2);
    }
  
  fprintf(stderr,"Loading parameters file...");
  fflush(stderr);
  
  if (fscanf(fp,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %d",
	     &FsP,&StartP,&StartB,&StopP,&AccP,&AccB,
	     &DonP,&DonB,&BlastS,&EstP,&EstM) != 11)
    {
      fprintf(stderr, "Incorrect parameter file EuGene.par\n");
      exit(2);
    }
  fprintf(stderr,"done\n");

  // any Frameshift prob below -1000.0 means "not possible"
  if (FsP <= -1000.0) FsP = NINFINITY;

  while ((carg = getopt(argc, argv, "dshbm:w:f:n:o:p:x:y:u:v:g:")) != -1) {
    
    switch (carg) {
            
    case 'n':            /* -n normalize across frames      */
      if (! GetIArg(optarg, &normopt, normopt))
	errflag++;
      break;
      
    case 'h':           /* help                             */
      errflag++;
      break;

    case 's':           /* Single gene mode. Sets Prio to 0/0/1   */
      ExonPrior = 0.0;
      IntronPrior = 0.0;
      InterPrior = 1.0;
      break;

    case 'c':            /* -c couverture      */
      if (! GetIArg(optarg, &golap, golap))
        errflag++;
      break;
      
    case 'l':            /* -l imglen      */
      if (! GetIArg(optarg, &glen, glen))
        errflag++;
      break;
      
    case 'u':            /* -u From      */
      if (! GetIArg(optarg, &gfrom, gfrom))
        errflag++;
      break;
      
    case 'v':            /* -v To      */
      if (! GetIArg(optarg, &gto, gto))
        errflag++;
      break;
      
    case 'x':            /* -x resx      */
      if (! GetIArg(optarg, &resx, resx))
        errflag++;
      break;
      
    case 'y':            /* -y resy      */
      if (! GetIArg(optarg, &resy, resy))
        errflag++;
      break;
      
    case 'g':           /* -g "Graphic File"                */
      strcpy(grname,optarg);
      graph = TRUE;
      break;
    
    case 'p':           /* print opt: short/long/detailed   */
      printopt = optarg[0];
      if ((printopt != 's') && (printopt != 'l') && (printopt != 'd'))
	errflag++;
      break;
      
    case 'm':           /* -m matrix                        */
      (void) strcpy(matname, optarg);
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
      blastopt = TRUE;
      break; 
      
    case 'd':           /* -d  use cDNA blastn results      */
      estopt = TRUE;
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
    fprintf(stderr, "usage: EuGene [-h] [-m matrix] [-b] [-c] [-p s|l|d] [-w window]\n");
    fprintf(stderr, "              [-n 0|1|2] [-o offset] [-s] [-f FS prob] FASTA files\n");
    exit(1);
  }
  
  // open matrix
  
  if (! (fp = OpenFile(matname,  "rb"))) {
    fprintf(stderr, "cannot open matrix file %s\n",  matname);
    exit(2);
  }

  fprintf(stderr,"Loading Markov model...");
  fflush(stderr);
  
  for  (i = 0;  i < 5;  i ++) {
    IMMatrix[i] = new BString_Array(MODEL_LEN, ALPHABET_SIZE);
    IMMatrix[i]->Read(fp);
    fprintf(stderr,"%d ",i+1);
    fflush(stderr);
  }
  fprintf(stderr,"done\n");

  int sequence;
  for (sequence = optind; sequence < argc ; sequence++) {
      
    (void) strcpy(fstname, argv[sequence]);
    
    // read fasta file
  
  Input_Size = INIT_SIZE;
  Data  = (char *) Safe_malloc(Input_Size);
  
  fp = (*fstname ? File_Open (fstname, "r") : stdin);
  
  if (fp == NULL) {
    fprintf(stderr, "cannot open fasta file %s\n",  fstname);
    exit(3);
  }

  fprintf(stderr,"Loading sequence...");
  fflush(stderr);
  
  Read_String (fp, Data, Input_Size, seqname, FALSE);

  Data += 1;
  Data_Len = strlen(Data);
  fprintf(stderr,"%s, %ld bases read\n",seqname, Data_Len);

  if (fp != stdin) fclose(fp);
  
  Ch_Ct[0] =  Ch_Ct [1] = Ch_Ct [2] = Ch_Ct [3] = 0.0;

  for  (i = 0;  i < Data_Len;  i ++)
    {
      Data [i] = Filter (tolower (Data [i]));
      switch  (Data [i])
        {
	case  'a' :
	case  't' :
	  Ch_Ct [0] += 1.0;
	  break;
	case  'c' :
	case  'g' :
	  Ch_Ct [1] += 1.0;
	  break;
        }
    }
  
  Ch_Ct [2] = Ch_Ct [1];
  Ch_Ct [3] = Ch_Ct [0];
  
  for  (i = 0;  i < 4;  i ++)
    Ch_Ct [i] = Ch_Ct [i] / (2.0 * Data_Len);
  
  fprintf (stderr,"GC Proportion = %.1f%%\n", 200.0 * Ch_Ct [1]);
  
  for  (i = 0;  i < 4;  i ++)
    Ch_Ct [i] = log (Ch_Ct [i]);

  for (i = 0;  i < 9;  i ++)
    BaseScore[i] = new REAL[Data_Len+1];
 
  Stop[0] = new unsigned char[Data_Len+1];
  Stop[1] = new unsigned char[Data_Len+1];

  ESTMatch = new unsigned char[Data_Len+1];

  for (i = 0; i<= Data_Len; i++) ESTMatch[i] = 0;

  Start[0] = new REAL[Data_Len+1];
  Start[1] = new REAL[Data_Len+1];

  Acc[0] = new REAL[Data_Len+1];
  Acc[1] = new REAL[Data_Len+1];

  Don[0] = new REAL[Data_Len+1];
  Don[1] = new REAL[Data_Len+1];
  
  Find_Stop_Codons(Data, Data_Len, Stop);

  fprintf(stderr, "Reading start files...");
  fflush(stderr);
  
  strcpy(tempname,fstname);
  strcat(tempname,".starts");
  Read_Start(tempname,Data_Len, exp(-StartP), StartB, Start[0],0);
  fprintf(stderr,"forward,");
  fflush(stderr);

  strcpy(tempname,fstname);
  strcat(tempname,".startsR");
  Read_Start(tempname,Data_Len, exp(-StartP), StartB, Start[1],1);
  fprintf(stderr," reverse done\n");

  fprintf(stderr, "Reading splice site files...");  
  fflush(stderr);

  Read_Splice(fstname,1,Data_Len, Acc[0], Don[0],AccP,AccB,DonP,DonB);
  fprintf(stderr,"forward,");
  fflush(stderr);

  Read_Splice(fstname,-1,Data_Len, Acc[1], Don[1],AccP,AccB,DonP,DonB);
  fprintf(stderr," reverse done\n");

  if (graph) {
    if (gfrom != -1) sprintf(grname+strlen(grname),"%d",gfrom);
    if (gto != -1) sprintf(grname+strlen(grname),"-%d",gto);
    
    if ((gfrom <= 0)|| (gfrom >= Data_Len)) gfrom = 1;
    if ((gto <= 0)  || (gto <= gfrom) || (gto > Data_Len)) gto = Data_Len;
    
    gfrom--;
    gto--;
    
    if (glen < 0) glen = ((gto-gfrom+1 <= 6000) ? gto-gfrom+1 : 6000);
    
    InitPNG(resx,resy,gfrom,gto,golap,glen,grname);
  }
    
  // Precompute all probabilities
  
  fprintf(stderr, "Computing coding probabilities...");  
  fflush(stderr);
  
  for  (i = 0 ; i < Data_Len ; i++)
    {
      Simple_Score (Data, i, Data_Len,
		    IMMatrix, Score);
            
      for (j = 0;  j < 9;  j ++)
	BaseScore[j][i] = Score[j];
    }

  for (j = 0;  j < 9;  j ++)
    BaseScore[j][Data_Len] = 1.0;

  fprintf(stderr,"done\n");



  // Blastx against protein databases 1st simple version, we simply
  // enhance coding probability according to phase Dangerous (NR
  // contains translation of frameshifted sequences) another
  // possibility would be too forbid introns/intergenic states
 
  if (blastopt)
    {
      FILE *fblast;
       int deb,fin,phase,score,Pfin,ProtDeb,ProtFin,PProtFin,PPhase;
       char A[128],B[128];
       char *ProtId, *PProtId,*tmp;
       REAL BScore;
       REAL *BlastScore[6];
       
       for (i=0; i<6; i++)
	 {
	   BlastScore[i] = new REAL[Data_Len];
	   for (j=0; j<Data_Len; j++)
	     BlastScore[i][j] = 0.0;
	 }
       
       strcpy(tempname,fstname);
       strcat(tempname,".blast");
       fblast = fopen(tempname, "r");
       
       if (fblast == NULL)
	 {
	   fprintf (stderr, "ERROR:  Could not open blast file \"%s\"\n", 
		    tempname);
	   exit (-1);
	 }

       
       A[0] = B[0]= 0;
       ProtId = A;
       PProtId = B;
       PProtFin = -10000;

       fprintf(stderr,"Reading Blastx data...");
       fflush(stderr);

       while (fscanf(fblast,
		     "%d %d %d %*s %d %s %d %d\n", 
		     &deb, &fin, &score, &phase,ProtId,&ProtDeb,&ProtFin) != EOF) 
	 {
	   if (phase < 0) {
	     j = deb;
	     deb = fin;
	     fin = j;
	     j = ProtDeb;
	     ProtDeb = ProtFin;
	     ProtFin = j;
	   }
	  
	   //	     printf("%s match (%d %d) - (%d %d) p. %d/%d, prev %d\n",ProtId,deb,fin,ProtDeb,ProtFin,PPhase,phase,PProtFin);
	   if ((strcmp(ProtId,PProtId) == 0) && (abs(PProtFin-ProtDeb) <= 8)) {
	     //	     printf("Sequential hit detected\n");

	     //	     for (i = Pfin; i< deb; i++)
	       //	       BaseScore[8][i] = (BaseScore[8][i]*100.0)/
	       //		 (100.0+BlastS);
	     if (graph) PlotLine(Pfin,deb,PPhase,phase,0.6,0.6,4);
	   }

	   Pfin = fin;
	   tmp = PProtId;
	   PProtId = ProtId;
	   ProtId = tmp;
	   PProtFin = ProtFin;
	   PPhase = phase;


	   phase = ph06(phase);
	   BScore = ((REAL)score)/((REAL)abs(fin-deb));
	   
	   for (i = deb-1; i < fin; i++)  
	     {
	       if (BScore > BlastScore[phase][i]) 
		 BlastScore[phase][i] = BScore;
	     }
	 }

       fprintf(stderr,"done\n");

       for (i=0; i<6; i++)
	 {
	   for (j=0; j<Data_Len; j++)
	     {
	       BaseScore[i][j] = 
		 ((BaseScore[i][j]*100.0)+(BlastScore[i][j]*BlastS))/
		 (100.0+(BlastScore[i][j]*BlastS));
	       BaseScore[i][12] = (BaseScore[i][12]*100.0)/
		 (100.0+(BlastScore[i][j]*BlastS));
	       if (graph && (BlastScore[i][j] > 0.0)) PlotBarI(j,PhaseAdapt(i),0.6,1,4);
		 
	     }
	   delete BlastScore[i];
	 }
       fclose(fblast);
    }
  

  /* Blastn against EST. 
     1st Simple version: we simply forbid introns during matches */

  if (estopt)
    {
      FILE *fblast;
      int deb,fin,Pfin,brin,EstDeb,EstFin,PEstFin;
      char A[128],B[128];
      char *EstId, *PEstId,*tmp;

       strcpy(tempname,fstname);
       strcat(tempname,".est");
       fblast = fopen(tempname, "r");
       
       if (fblast == NULL)
	 {
	   fprintf (stderr, "ERROR:  Could not open EST blast file \"%s\"\n", 
		    tempname);
	   exit (-1);
	 }
       
       A[0] = B[0]= 0;
       EstId = A;
       PEstId = B;
       PEstFin = -10000;

       fprintf(stderr,"Reading cDNA hits...");
       fflush(stderr);

       while (fscanf(fblast,
		     "%d %d %*s %*s %d %s %d %d\n", 
		     &deb, &fin,&brin,EstId,&EstDeb,&EstFin) != EOF)
	 {
	   //	   printf("seq %s  matche (%d %d) sur (%d %d), prev %d\n",EstId,deb,fin,EstDeb,EstFin,PEstFin);

	   // Use 1 or -1 to encode EST hit. Sign indicates forward/reverse
	   for (i = deb-1+EstM; i < fin-EstM; i++)  ESTMatch[i] |= (HitForward << (brin*4));
	   //	     if (!(ESTMatch[i] & NotAHit)) 

	   if ((strcmp(EstId,PEstId) == 0) && (abs(PEstFin-EstDeb) <= 6)) {

	     //	     printf("Sequential hit detected\n");
	     if (deb > 2*EstM) {

	     for (i = Pfin-EstM; i < Pfin-1+EstM; i++) 
	       ESTMatch[i] = (ESTMatch[i] & (0xf0 >> (brin*4))) | (MLeftForward << (brin*4));

	     for (i = Pfin-1+EstM; i < deb-EstM; i++)  
	       ESTMatch[i] = (ESTMatch[i] & (0xf0 >> (brin*4))) | (GapForward << (brin*4));

	     for (i = deb-EstM; i < deb-1+EstM; i++)  
	       ESTMatch[i] = (ESTMatch[i] & (0xf0 >> (brin*4))) | (MRightForward << (brin*4));

	     if (graph) PlotLine(Pfin-1+EstM,deb-EstM,4-(brin*8),4-(brin*8),0.8,0.8,2);
	     }
	   }
	   
	   Pfin = fin;
	   tmp = PEstId;
	   PEstId = EstId;
	   EstId = tmp;
	   PEstFin = EstFin;
	 }
	     
       deb = 0;
       for (i=0; i<Data_Len; i++) {
	 if (graph && (ESTMatch[i] & HitForward)) PlotBarI(i, 4,0.6,1,2);
	 if (graph && (ESTMatch[i] & HitReverse)) PlotBarI(i, -4,0.6,1,2);
	 if (Inconsistent(ESTMatch[i])) {
	   deb = 1;
	 }
       }
       fclose(fblast);

       fprintf(stderr,"done\n");
       if (deb) 
	 fprintf(stderr,"  -> Warning: inconsistent cDNA data !\n");
    }

  // Data allocation for the shortest path with length constraints algorithm
  //
  // 0  codant 1
  // 1  codant 2
  // 2  codant 3
  // 3  codant -1
  // 4  codant -2
  // 5  codant -3
  // 6  Intron frame 0
  // 7  Intron frame 1
  // 8  Intron frame 2
  // 9  Intron frame 0
  // 10 Intron frame 1
  // 11 Intron frame 2  
  // 12 intergenique
  
  char *Choice;
  BackPoint *LBP[13];
  double BestU;
  char best;
  unsigned char Switch;

  Choice =  new char[Data_Len+2];
  
  for (i=0; i<13; i++) {
    LBP[i] = new BackPoint(i, -1,0.0);
    LBP[i]->Next = LBP[i];
    LBP[i]->Prev = LBP[i];
  }
  
  // prior on the initial state  
  // Selon Sato et al 1999 / Terryn et al. 1999 
  
  LBP[0]->Update(log(ExonPrior/6.0)/2.0);
  LBP[1]->Update(log(ExonPrior/6.0)/2.0);
  LBP[2]->Update(log(ExonPrior/6.0)/2.0);
  LBP[3]->Update(log(ExonPrior/6.0)/2.0);
  LBP[4]->Update(log(ExonPrior/6.0)/2.0);
  LBP[5]->Update(log(ExonPrior/6.0)/2.0);
  
  LBP[6]->Update(log(IntronPrior/6.0)/2.0);
  LBP[7]->Update(log(IntronPrior/6.0)/2.0);
  LBP[8]->Update(log(IntronPrior/6.0)/2.0);
  LBP[9]->Update(log(IntronPrior/6.0)/2.0);
  LBP[10]->Update(log(IntronPrior/6.0)/2.0);
  LBP[11]->Update(log(IntronPrior/6.0)/2.0);
  
  LBP[12]->Update(log(InterPrior)/2.0);
  
  REAL  maxi, PBest[23];
  BackPoint *PrevBP[23];

  for (i = 0; i <= Data_Len; i++) {
    
    for (k = 0 ; k < 19; k++) 
      PrevBP[k] = LBP[k%13]->BestUsable(i,SwitchMask[k],MinLength[k],&PBest[k]);

    for (k = 0; k < 4; k++) 
      PrevBP[k+19] = LBP[12]->StrictBestUsable(i,StrictSwitchMask[k],StrictMinLength[k],&PBest[k+19]);

    //    for (k = 0; k < 6; k++) 
    //      if (PBest[k] != PBest[k+13]) printf("pos %d state %d %f <> %f\n",i,k,PBest[k],PBest[k+13]);
//     BestU = MAX(PBest[19],PBest[20]);
//     if (BestU != PBest[12]) printf("pos %d state %d %f <> MAX(%f,%f)\n",
// 				   i,k,PBest[12],PBest[19],PBest[20]);
//     BestU = MAX(PBest[21],PBest[22]);
//     if (BestU != PBest[12]) printf("pos %d state %d %f <> MAX(%f,%f)\n",
// 				   i,k,PBest[12],PBest[21],PBest[22]);

    // Codant en forward
    for (k = 0; k<3; k++) {
      
      maxi = NINFINITY;     

      // S'il y  a un STOP en phase on ne peut continuer
      if ((i % 3 == k) && Stop[0][i]) 
	LBP[k]->Update(NINFINITY); 
      else 
	LBP[k]->Update(log(1.0-Don[0][i]));
      
      LBP[k]->BestUsable(i,SwitchAny,0,&BestU);
      if (BestU > maxi) {
	maxi = BestU;
	best = k;
      }
      
      // - on commence a coder (Start)

      // Il y a deux possibilites suivant le signal qui est utilise
      // pour demarrer la zone intergenique:
      //
      //   - un STOP : on a deux genes confluents. La contrainte de
      //     distance est MinIgFlow -  cas 19
      //
      //   - un START : on a deux genes divergents. La contrainte de
      //     distance est MinIgDiv. - cas 20
      //
      // On prend le meilleur des deux.

      if ((i % 3 == k) && Start[0][i] != 0.0) {
	BestU = PBest[19]+log(Start[0][i]);
	if (BestU > maxi) {
	  maxi = BestU;
	  best = 19;
	  Switch = SwitchStart;
	}
	BestU = PBest[20]+log(Start[0][i]);
	if (BestU > maxi) {
	  maxi = BestU;
	  best = 20;
	  Switch = SwitchStart;
	}
      }
      
      // - on recommence a coder (Accepteur)
      BestU = PBest[6+((i-k+3) % 3)]+log(Acc[0][i]);
      if (BestU > maxi) {
	maxi = BestU;
	best = 6+((i-k+3) % 3);
	Switch = SwitchAcc;
      }
      
      // - il y a insertion
      //      b[3] = LR[i%2][(k+2)%3]+FsP;
      //      a[3] = (k+2)%3;
      
      // - il y a deletion
      //      b[4] = LR[i%2][(k+1)%3]+FsP;
      //      a[4] = (k+1)%3;
      

      if (best != k) LBP[k]->InsertNew(((best >= 19) ? 12 : best),Switch,i,maxi,PrevBP[best]);
      
      LBP[k]->Update(log(BaseScore[k][i])+((ESTMatch[i] & Gap) != 0)*EstP);
    }
      
    // Codant en reverse
    
    for (k = 3; k<6; k++) {

      maxi = NINFINITY;

      // On continue sauf si l'on rencontre un autre STOP
      if (((Data_Len-i) % 3 == k-3) && Stop[1][i]) 
	LBP[k]->Update(NINFINITY);
      else 
	LBP[k]->Update(log(1.0-Don[1][i]));
      
      LBP[k]->BestUsable(i,SwitchAny,0,&BestU);
      if (BestU > maxi) {
	maxi = BestU;
	best = k;
      }

      // - on commence a coder (Stop)
      
      // Il y a deux possibilites suivant le signal qui est utilise
      // pour demarrer la zone intergenique:
      //
      //   - un START : on a deux genes confluents. La contrainte de
      //     distance est MinIgFlow -  cas 21
      //
      //   - un STOP : on a deux genes convergents. La contrainte de
      //     distance est MinIgConv -  cas 22
      //
      // On prend le meilleur des deux.

      if (((Data_Len-i) % 3 == k-3) && Stop[1][i]) {// && !(ESTMatch[i] & Gap)) {
	BestU = PBest[21]-StopP;
	if (BestU > maxi) {
	  maxi = BestU;
	  best = 21;
	  Switch = SwitchStop;
	}
	BestU = PBest[22]-StopP;
	if (BestU > maxi) {
	  maxi = BestU;
	  best = 22;
	  Switch = SwitchStop;
	}
      }

      // - on recommence a coder (Donneur)
      BestU = PBest[9+((Data_Len-i-k) % 3)]+log(Don[1][i]);
      if (BestU > maxi) {
	maxi = BestU;
	best = 9+((Data_Len-i-k) % 3);
	Switch =  SwitchDon;
      }
      
      // - il y a insertion
      ///	b[3] = LR[i%2][((k+2) % 3)+3]+FsP;
      //	a[3] =((k+2) % 3)+3 ;
      
      // - il y a deletion
      //	b[4] = LR[i%2][((k+1) % 3)+3]+FsP;
      //	a[4] = ((k+1) % 3)+3;
      
      if (best != k) LBP[k]->InsertNew(((best >= 19) ? 12 : best),Switch,i,maxi,PrevBP[best]);
      
      LBP[k]->Update(log(BaseScore[k][i]) + ((ESTMatch[i] & Gap) != 0)*EstP);
    }

    // Intergenique

    maxi = NINFINITY;

    // On reste intergenique
    LBP[12]->Update(log(1.0-Start[0][i])+log(1.0-Start[1][i]));
    
    LBP[12]->BestUsable(i,SwitchAny,0,&BestU);
    if (BestU > maxi) {
      maxi = BestU;
      best = 12;
    }
 
    // Forward
    for (k = 0; k<3; k++) {
      // - on quitte un codant
      // y a-t-il un Stop en phase ?
      if ((i % 3 == k) && Stop[0][i]) {// && !(ESTMatch[i] & Gap)) {
	BestU = PBest[k+13]-StopP;
	if (BestU > maxi) {
	  maxi = BestU;
	  best = k+13;
	  Switch = SwitchStop;
	}
      }
    }

    // Reverse
    for (k = 3; k < 6; k++) {
      // - on quitte un codant
      // y a-t-il un Start en phase ?
      if (((Data_Len-i) % 3 == k-3) && Start[1][i] != 0.0) {
	BestU = PBest[k+13]+log(Start[1][i]);
	if (BestU > maxi) {
	  maxi = BestU;
	  best = k+13;
	  Switch = SwitchStart;
	}
      }
    }
      
    if (best != 12) LBP[12]->InsertNew(best % 13,Switch,i,maxi,PrevBP[best]);
      
    LBP[12]->Update(log(BaseScore[8][i]));
      
    // Introns de phase k (Forward)
    for (k = 0; k<3; k++) {

      maxi = NINFINITY;

      // On reste intronique
      LBP[6+k]->Update(log(1.0-Acc[0][i]));
      LBP[6+k]->BestUsable(i,SwitchAny,0,&BestU);
      if (BestU > maxi) {
	maxi = BestU;
	best = 6+k;
      }
      // - on quitte un exon
      BestU = PBest[((i-k+3) % 3)]+log(Don[0][i]);
      if (BestU > maxi) {
	maxi = BestU;
	best = ((i-k+3) % 3);
	Switch = SwitchDon;
      }
 
      if (best != 6+k) LBP[6+k]->InsertNew(best,Switch,i,maxi,PrevBP[best]);
      
      LBP[6+k]->Update(log(BaseScore[6][i])+((ESTMatch[i] & Hit) != 0)*EstP);
    }
      
    // Introns de phase -k
    for (k = 0; k<3; k++) {

      maxi = NINFINITY;

      // On reste intronique
      LBP[9+k]->Update(log(1.0-Acc[1][i]));
      LBP[9+k]->BestUsable(i,SwitchAny,0,&BestU);
      if (BestU > maxi) {
	maxi = BestU;
	best = 9+k;
      }

      // - on quitte un exon
      BestU = PBest[3+((Data_Len-i-k) % 3)]+log(Acc[1][i]);
      if (BestU > maxi) {
	maxi = BestU;
	best = 3+((Data_Len-i-k) % 3);
	Switch = SwitchAcc;
      }
      
      
      if (best != 9+k) LBP[9+k]->InsertNew(best,Switch,i,maxi,PrevBP[best]);
      
      LBP[9+k]->Update(log(BaseScore[7][i])+((ESTMatch[i] & Hit) != 0)*EstP);
    }
  }

  LBP[0]->Update(log(ExonPrior/6.0)/2.0);
  LBP[1]->Update(log(ExonPrior/6.0)/2.0);
  LBP[2]->Update(log(ExonPrior/6.0)/2.0);
  LBP[3]->Update(log(ExonPrior/6.0)/2.0);
  LBP[4]->Update(log(ExonPrior/6.0)/2.0);
  LBP[5]->Update(log(ExonPrior/6.0)/2.0);
  LBP[6]->Update(log(IntronPrior/6.0)/2.0);
  LBP[7]->Update(log(IntronPrior/6.0)/2.0);
  LBP[8]->Update(log(IntronPrior/6.0)/2.0);
  LBP[9]->Update(log(IntronPrior/6.0)/2.0);
  LBP[10]->Update(log(IntronPrior/6.0)/2.0);
  LBP[11]->Update(log(IntronPrior/6.0)/2.0);
  LBP[12]->Update(log(InterPrior)/2.0);

  for (i = 0; i < 13; i++) {
    PrevBP[i] = LBP[i]->BestUsable(Data_Len+1,SwitchAny,0,&BestU);    
    LBP[i]->InsertNew(i,0xFF,Data_Len+1,BestU,PrevBP[i]);
  }
  
  // Backtrace in the graph
  LBP[0]->BestUsable(Data_Len+1,SwitchAny,0,&maxi);
  j = 0;

  for (i = 1; i < 13 ; i++) {
    LBP[i]->BestUsable(Data_Len+1,SwitchAny,0,&BestU);
    if (BestU > maxi) {
      maxi = BestU;
      j = i;
      }
  }

  LBP[j]->BackTrace(Choice);
  
  for  (i = 0;  i < 13;  i ++) LBP[i]->Zap();


  fprintf(stderr,"Optimal path length = %#f\n",maxi);

  if (graph) PlotPredictions(Choice,BaseScore,Stop,Start,Acc,Don);

#include "Output.h"

  // free used memory
  if (graph) ClosePNG();

  delete (Data-1);
    
  delete Stop[0];
  delete Stop[1];

  delete ESTMatch;
  
  delete Start[0];
  delete Start[1];

  delete Don[0];
  delete Don[1];
  
  delete Acc[0];
  delete Acc[1];

  delete Choice;

  for  (i = 0;  i < 9;  i ++) delete BaseScore[i];
  }

  // free remaining used memory
  
  for  (i = 0;  i < 5;  i ++)
    delete IMMatrix[i];
    
  return  0;
}
