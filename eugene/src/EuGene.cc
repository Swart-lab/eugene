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
// ------------------ Globals ---------------------
const int  MODEL_LEN = 9;
const int  SIMPLE_MODEL_LEN = 6;
const int  ALPHABET_SIZE = 4;
const REAL NINFINITY = log(0.0);
const int PredWidth = 2;
const REAL Epsilon = 1e-2;

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
#ifdef UseSpliceP
#include "SpliceP.h"
#else
#include "SpliceN.h"
#endif
#include "../GDIF/gdIF.h"

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
  int i,j,p,prev = 0;

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
	PlotBarF(i,4,0.5,NORM(log(Acc[0][i]),16),4);
      
      if (Acc[1][i+1] > 0.0) 
	PlotBarF(i,-4,0.5,NORM(log(Acc[1][i+1]),16),4);

      if (Don[0][i] > 0.0) 
	PlotBarF(i,4,0.5,NORM(log(Don[0][i]),16),5);
      
      if (Don[1][i+1] > 0.0) 
	PlotBarF(i,-4,0.5,NORM(log(Don[1][i+1]),16),5);
    }    
}

// -------------------------------------------------------------------------
//            MAIN
// -------------------------------------------------------------------------
main  (int argc, char * argv [])
{
  int              i, j, k, carg, errflag;
  FILE             *fp;
  REAL             *BaseScore[9],*Start[2] ;
  REAL             *Don[2],*Acc[2];
  unsigned char    *Stop[2],*NoStop;
  long int         Input_Size;
  REAL             Ch_Ct[ALPHABET_SIZE] = {0.0};
  BString_Array    *IMMatrix[5];
  char             *Data;
  char             printopt;
  char             matname[FILENAME_MAX+1], tempname[FILENAME_MAX+1], fstname[FILENAME_MAX+1];
  char             *grname;
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
      grname = optarg;
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
      
    case 'd':           /* -c  use cDNA blastn results      */
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
  for (sequence = optind; sequence < argc ; sequence++)
    {

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
  fprintf(stderr,"%d bases read\n",Data_Len);

  if (gfrom < 0) gfrom = 0;
  if ((gto == -1) || (gto >= Data_Len))   gto = Data_Len-1;
  if (glen < 0) glen = ((Data_Len < 6000) ? Data_Len : 6000);

  if (fp != stdin)
    fclose(fp);
  
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
    BaseScore[i] = (REAL *)Safe_malloc(sizeof(REAL)*(Data_Len+1));
 
  Stop[0] = (unsigned char *)Safe_malloc(sizeof(unsigned char)*(Data_Len+1));
  Stop[1] = (unsigned char *)Safe_malloc(sizeof(unsigned char)*(Data_Len+1));

  NoStop = (unsigned char *)Safe_malloc(sizeof(unsigned char)*(Data_Len+1));

  for (i = 0; i<= Data_Len; i++) NoStop[i] = 0;

  Start[0] = (REAL *)Safe_malloc(sizeof(REAL)*(Data_Len+1));
  Start[1] = (REAL *)Safe_malloc(sizeof(REAL)*(Data_Len+1));

  Acc[0] = (REAL *)Safe_malloc(sizeof(REAL)*(Data_Len+1));
  Acc[1] = (REAL *)Safe_malloc(sizeof(REAL)*(Data_Len+1));

  Don[0] = (REAL *)Safe_malloc(sizeof(REAL)*(Data_Len+1));
  Don[1] = (REAL *)Safe_malloc(sizeof(REAL)*(Data_Len+1));
  
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

  if (graph) InitPNG(resx,resy,gfrom,gto,golap,glen,grname);

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
	   BlastScore[i] = (REAL *)Safe_malloc(sizeof(REAL)*(Data_Len));
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
       
       while (fscanf(fblast,
		     "%d %d %d %*s %d %s %d %d\n", 
		     &deb, &fin, &score, &phase,ProtId,&ProtDeb,&ProtFin) != EOF) 
	 {
	   if (phase < 0) {
	     j = deb;
	     deb = fin;
	     fin = j;
	   }
	  
	   if ((strcmp(ProtId,PProtId) == 0) && (abs(PProtFin-ProtDeb) <= 6) && graph) {
	     PlotLine(Pfin,deb,PPhase,phase,0.6,0.6,4);
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
       
       for (i=0; i<6; i++)
	 {
	   for (j=0; j<Data_Len; j++)
	     {
	       BaseScore[i][j] = 
		 ((BaseScore[i][j]*100.0)+(BlastScore[i][j]*BlastS))/
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
      char *ESTHit[2];
      int deb,fin,Pfin,brin,EstDeb,EstFin,PEstFin;
      char A[128],B[128];
      char *EstId, *PEstId,*tmp;

      for (i=0; i<2; i++)
	{
	  ESTHit[i] = (char *)Safe_malloc(sizeof(char)*(Data_Len));
	  for (j = 0; j < Data_Len; j++)
	    ESTHit[i][j] = 0;
	}
      
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

       while (fscanf(fblast,
		     "%d %d %*s %*s %d %s %d %d\n", 
		     &deb, &fin,&brin,EstId,&EstDeb,&EstFin) != EOF)
	 {
	   for (i = deb-1+EstM; i < fin-EstM; i++)  ESTHit[brin][i] = 1;
	   if ((strcmp(EstId,PEstId) == 0) && (abs(PEstFin-EstDeb) <= 6)) {
	     for (i = Pfin-1+EstM; i < deb-EstM; i++)  NoStop[i] = 1;
	     if (graph) PlotLine(Pfin-1+EstM,deb-EstM,4-(brin*8),4-(brin*8),0.6,0.6,2);
	   }
	   Pfin = fin;
	   tmp = PEstId;
	   PEstId = EstId;
	   EstId = tmp;
	   PEstFin = EstFin;
	 }
	     

       for (i=0; i<Data_Len; i++)
	 {
	   if (ESTHit[0][i]) {
	     BaseScore[6][i] = Epsilon;
	     if (graph) PlotBarI(i, 4,0.6,1,2);
	   }
	   if (ESTHit[1][i]) {
	     BaseScore[7][i] = Epsilon;
	     if (graph) PlotBarI(i,-4,0.6,1,2);
	   }
	 }
       delete ESTHit[0];
       delete ESTHit[1];
       fclose(fblast);
    }

  // Data allocation for Belmann
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
  
  char *macs[13], *Choice;
  REAL *LR[2];

  LR[0] = (REAL *)Safe_malloc(sizeof(REAL)*13);
  LR[1] = (REAL *)Safe_malloc(sizeof(REAL)*13);
  Choice =  (char *)Safe_malloc(sizeof(char)*(Data_Len+1));

  for (i=0; i<13; i++)
    {
      macs[i] = (char *)Safe_malloc(sizeof(char)*(Data_Len+2));
      macs[i][0] = -1;
      LR[0][i] = LR[1][i] = 0.0;
    }

  // prior on the initial state  
  // Selon Sato et al 1999 / Terryn et al. 1999 



  LR[0][0] =  log(ExonPrior/6.0)/2.0;
  LR[0][1] =  log(ExonPrior/6.0)/2.0;
  LR[0][2] =  log(ExonPrior/6.0)/2.0;
  LR[0][3] =  log(ExonPrior/6.0)/2.0;
  LR[0][4] =  log(ExonPrior/6.0)/2.0;
  LR[0][5] =  log(ExonPrior/6.0)/2.0;
  
  LR[0][6] =  log(IntronPrior/6.0)/2.0;
  LR[0][7] =  log(IntronPrior/6.0)/2.0;
  LR[0][8] =  log(IntronPrior/6.0)/2.0;
  LR[0][9] =  log(IntronPrior/6.0)/2.0;
  LR[0][10] = log(IntronPrior/6.0)/2.0;
  LR[0][11] = log(IntronPrior/6.0)/2.0;

  LR[0][12] = log(InterPrior)/2.0;

  REAL  b[7],m;
  char    a[7];

  // on ne compte plus les BaseScore a chaque fois car cela ne peut
  // affecter l'ordre. On l'ajoute simplement apres coup a LR

  for (i = 0; i <= Data_Len; i++)
    {
      // Codant en forward
      for (k = 0; k<3; k++)
	{
	  // S'il y  a un STOP en phase c'est fini
	  if ((i % 3 == k) && Stop[0][i])
	    b[0] = NINFINITY;
	  else b[0] = LR[i%2][k]+log(1.0-Don[0][i]);
	  a[0] = k;
	  
	  // - on commence a coder (Start)
	  b[1] = LR[i%2][12]+
	    ((i % 3 == k) ? log(Start[0][i]) : NINFINITY);
	  a[1] = 12;
	  
	  // - on recommence a coder (Accepteur)
	  b[2] = LR[i%2][6+((i-k+3) % 3)]+
	    ((i > 0) ? log(Acc[0][i-1]) : NINFINITY);
	  a[2] = 6+((i-k+3) % 3);

	  // - il y a insertion
	  b[3] = LR[i%2][(k+2)%3]+FsP;
	  a[3] = (k+2)%3;
	  
	  // - il y a deletion
	  b[4] = LR[i%2][(k+1)%3]+FsP;
	  a[4] = (k+1)%3;
	  
	  m = b[0];
	  macs[k][i+1] = a[0];
	  
	  for (j = 1; j < 5; j++)
	    if (b[j] > m)
	      {
		m = b[j];
		macs[k][i+1] = a[j];
	      }
	  LR[(i+1)%2][k] = m+log(BaseScore[k][i]);
	}

      // Codant en reverse

      for (k = 3; k<6; k++)
	{
	  // On continue sauf si l'on rencontre un autre STOP
	  if (((Data_Len-i) % 3 == k-3) && Stop[1][i])
	    b[0] = NINFINITY;
	  else
	    b[0] = LR[i%2][k]+((i < Data_Len) ? log(1.0-Don[1][i]) : 0.0);

	  a[0] = k;
	  
	  // - on commence a coder (Stop)
	  b[1] = LR[i%2][12]+
	    ((((Data_Len-i) % 3 == k-3) && Stop[1][i] && !NoStop[i]) ?
	     (-StopP) : NINFINITY);
	  a[1] = 12;
	  
	  // - on recommence a coder (Donneur)
	  b[2] = log(Don[1][i])+LR[i%2][9+((Data_Len-i-k) % 3)];
	  a[2] = 9+((Data_Len-i-k) % 3);

          // - il y a insertion
         b[3] = LR[i%2][((k+2) % 3)+3]+FsP;
         a[3] =((k+2) % 3)+3 ;
 
         // - il y a deletion
         b[4] = LR[i%2][((k+1) % 3)+3]+FsP;
         a[4] = ((k+1) % 3)+3;
	  
	  m = b[0];
	  macs[k][i+1] = a[0];

	  for (j = 1; j < 5; j++)
	    if (b[j] > m)
	      {
		m = b[j];
		macs[k][i+1] = a[j];
	      }
	  LR[(i+1)%2][k] = m+log(BaseScore[k][i]);
	}

      // Intergenique
      // Forward
      for (k = 0; k<3; k++)
	{
	  // - on quitte un codant
	  b[k] = LR[i%2][k]+
	    // y a-t-il un Stop en phase ?
	    (((i % 3 == k) && Stop[0][i] && !NoStop[i]) ?
	     (-StopP) : NINFINITY);
	  a[k] = k;
	}

      // Reverse
      for (k = 3; k<6; k++)
	{
	  // - on quitte un codant
	  b[k] = LR[i%2][k]+
	    // y a-t-il un Start en phase ?
	    (((Data_Len-i) % 3 == k-3) ?
	     log(Start[1][i]) : NINFINITY);
	  a[k] = k;
	}
      
      // On reste intergenique
      b[6] = LR[i%2][12]+log(1.0-Start[0][i])+log(1.0-Start[1][i]);
      a[6] = 12;
      
      m = b[0];
      macs[12][i+1] = a[0];

      for (j = 1; j < 7; j++)
	if (b[j] > m)
	  {
	    m = b[j];
	    macs[12][i+1] = a[j];
	  }

      LR[(i+1)%2][12] = m+log(BaseScore[8][i]);


      // Introns de phase k (Forward)
      for (k = 0; k<3; k++)
	{
	  // - on quitte un exon
	  b[0] = LR[i%2][((i-k+3) % 3)]+log(Don[0][i]);
	  a[0] = ((i-k+3) % 3);
	  
	  // On reste intronique
	  b[1] = LR[i%2][6+k]+((i > 0) ? log(1.0-Acc[0][i-1]) : 0.0);
	  a[1] = 6+k;
	  
	  m = b[0];
	  macs[6+k][i+1] = a[0];
	  
	  if (b[1] > m)
	    {
	      m = b[1];
	      macs[6+k][i+1] = a[1];
	    }
	  LR[(i+1)%2][6+k] = m+log(BaseScore[6][i]);
	}
      
      // Introns de phase -k
      for (k = 0; k<3; k++)
	{
	  // - on quitte un exon
	  b[0] = LR[i%2][3+((Data_Len-i-k) % 3)]+
	    ((i < Data_Len) ? log(Acc[1][i+1]) : NINFINITY);
	  a[0] = 3+((Data_Len-i-k) % 3);
      
	  // On reste intronique
	  b[1] = LR[i%2][9+k]+log(1.0-Acc[1][i+1]);
	  a[1] = 9+k;

	  m = b[0];
	  macs[9+k][i+1] = a[0];

	  if (b[1] > m)
	    {
	      m = b[1];
	      macs[9+k][i+1] = a[1];
	    }
	  LR[(i+1)%2][9+k] = m+log(BaseScore[7][i]);
	}
    }

  LR[(Data_Len+1)%2][0]  +=  log(ExonPrior/6.0)/2.0;
  LR[(Data_Len+1)%2][1]  +=  log(ExonPrior/6.0)/2.0;
  LR[(Data_Len+1)%2][2]  +=  log(ExonPrior/6.0)/2.0;
  LR[(Data_Len+1)%2][3]  +=  log(ExonPrior/6.0)/2.0;
  LR[(Data_Len+1)%2][4]  +=  log(ExonPrior/6.0)/2.0;
  LR[(Data_Len+1)%2][5]  +=  log(ExonPrior/6.0)/2.0;
  
  LR[(Data_Len+1)%2][6]  +=  log(IntronPrior/6.0)/2.0;
  LR[(Data_Len+1)%2][7]  +=  log(IntronPrior/6.0)/2.0;
  LR[(Data_Len+1)%2][8]  +=  log(IntronPrior/6.0)/2.0;
  LR[(Data_Len+1)%2][9]  +=  log(IntronPrior/6.0)/2.0;
  LR[(Data_Len+1)%2][10] +=  log(IntronPrior/6.0)/2.0;
  LR[(Data_Len+1)%2][11] +=  log(IntronPrior/6.0)/2.0;

  LR[(Data_Len+1)%2][12] +=  log(InterPrior)/2.0;

  // Backtrace in the graph
  m = LR[(Data_Len+1)%2][0];
  j = 0;

  for (i = 1; i<13 ; i++)
    if (LR[(Data_Len+1)%2][i] > m)
      {
	m = LR[(Data_Len+1)%2][i];
	j = i;
      }
  
  for (i = Data_Len+1; i > 0; i--)
    {
      Choice[i-1] = j;
      j = macs[j][i];
    }
  
  fprintf(stderr,"Optimal path length = %#lf\n",m);

  if (graph) 
    {
      PlotPredictions(Choice,BaseScore,Stop,Start,Acc,Don);
    }

#include "Output.h"

  // free used memory
  if (graph) ClosePNG();

  for (i = 0; i< 13; i++)
    delete macs[i]; 

  delete (Data-1);
    
  delete Stop[0];
  delete Stop[1];
  
  delete Start[0];
  delete Start[1];

  delete Don[0];
  delete Don[1];
  
  delete Acc[0];
  delete Acc[1];

  delete Choice;
  delete LR[0];
  delete LR[1];

  for  (i = 0;  i < 9;  i ++)
    delete BaseScore[i];
    }

  // free remaining used memory
  
  for  (i = 0;  i < 5;  i ++)
    delete IMMatrix[i];
    
  return  0;
}
