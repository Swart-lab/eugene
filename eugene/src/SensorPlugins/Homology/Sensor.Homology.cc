#include "Sensor.Homology.h"

/*************************************************************
 **                      SensorHomology                     **
 *************************************************************/
// ---------------------------------------------------------------------------
// Travail sur la prise en compte des informations d'homologies.
// A partir de hits tblastx, on tente d'aider la prediction.
// Nouveau parametre TblastxP dans le fichier .par
// FORMATS necessaire:
// fichier matrice proteique: fichier BLOSUM ou PAM classique
// fichier .tblastx : comme les .blast, mais derniere colonne= seq prot. du hit subject.
// Ex:
// 802 864 104 0.0 +1 AB015498_354_416 354 416 PQGQTPLFPRIFGHEAGG*EKVGLWRVLEKV*LILHQEI

extern Parameters PAR;

// ----------------------
// fonction temporaire
// etude de la prise en compte possible des tblastx
// ----------------------
REAL SensorHomology :: tblastxupdate (int hitnb, REAL hitscore, double pen, double base) {
//  return (pen); // test avec score constant pour chq nt d'un hit
	// test d'integration de la donnee nbre de hits:
//	return (hitnb);    
//	return (hitnb+pen);
//	return (hitnb*base +pen);
	// test d'integration de la donnee score du hit:
//	return (hitscore);
//	return (hitscore +pen);
	return (hitscore*base +pen);
}

// ----------------------
//  Default constructor.
// ----------------------
SensorHomology :: SensorHomology(int n) : Sensor(n)
{
  TblastxP= PAR.getD("Homology.TblastxP");
  TblastxB= PAR.getD("Homology.TblastxB");
}

// ----------------------
//  Default destructor.
// ----------------------
SensorHomology :: ~SensorHomology ()
{
  for(int i=0;i<6;i++){
    delete TblastxNumber[i];
    delete TblastxScore[i];
  }
  delete [] TblastxNumber;
  delete [] TblastxScore;
}
// ----------------------
//  Init Homol.
// ----------------------
void SensorHomology :: Init (DNASeq *X)
{
  int i, j;
  int Len = X->SeqLen;
  
  type = Type_Content;
  
  FILE *ftblastx;
  FILE *protmatfile;
  int deb,fin,phase,ProtDeb,ProtFin,sens,score=0;
  float bits;
  REAL GlobalBits;
  const int MaxHitLen  = 15000;
  char tampon;
  char* paire= new char[3];
  paire[0]=paire[1]= '0';
  paire[2]='\0';
  ProtMat* PROTMAT;
  char tempname[FILENAME_MAX+1];
  
  TblastxNumber = new int*[6];
  TblastxScore = new REAL*[6];
  for (j = 0; j < 6; j++) {
    TblastxNumber[j]= new int[Len+1];
    TblastxScore[j] = new REAL[Len+1];
    for (i = 0; i<= Len; i++) {
      TblastxNumber[j][i]=0;
      TblastxScore[j][i]=0.0;
    }
  }
  
  strcpy(tempname,PAR.getC("Homology.protmatname"));
  // Lecture de la matrice proteique ("BLOSUM62" par defaut)
  protmatfile=FileOpen(NULL, tempname, "rt");
  if (protmatfile == NULL) {
    fprintf (stderr,"\ncannot open protein matrix file %s\n",tempname); exit (2);
  }
  // chargement classe protmat
  fprintf(stderr,"Loading protmat file.........");
  fichier2protmat(protmatfile , PROTMAT);
  fclose(protmatfile);
  //   PROTMAT->affichage(2);
  fprintf(stderr,"done\n");
  
  fprintf(stderr,"Reading tblastx data.........");
  fflush(stderr);
  strcpy(tempname,PAR.getC("fstname"));
  strcat(tempname,".tblastx");
  ftblastx = FileOpen(NULL,tempname, "r");
  if (ftblastx == NULL) {
    fprintf (stderr,"\ncannot open tblastx file %s\n",tempname); exit (2);
  }
  
  while (!feof(ftblastx)) {
    if (fscanf(ftblastx,"%d %d %f %*s %d %*s %d %d",
	       &deb, &fin, &bits, &phase, &ProtDeb,&ProtFin) != 6){
      if (feof(ftblastx)) break;
      fprintf(stderr,"error in tblastxfile format %s\n",tempname);
      exit(2);
    }
    if (abs(fin-deb) > MaxHitLen) {
      fprintf(stderr,"Similarity of extreme length rejected. Check tblastx file %s\n",tempname);
      continue;
    }
    sens=1;
    if (phase < 0) {
      sens = -1;
      j = deb;
      deb = fin;
      fin = j;
      j = ProtDeb;
      ProtDeb = ProtFin;
      ProtFin = j;
    }
    phase = ph06(phase);
    
    tampon=fgetc(ftblastx); // lecture de la suite (sequence hit subject)
    while (isspace(tampon)) tampon=fgetc(ftblastx); // en cas de plsrs separateurs avant la derniere colonne
    for (i = deb-1; i < fin; i++)  {
      if ( (feof(ftblastx)) ||
	   ( (isspace(tampon)) && (i != (fin)) ) ) {
	fprintf(stderr,"error in tblastx format, file %s\n",tempname);
	exit(2);
      }
      if ( (i-deb+1)%3 == 0 ) {
	if (i >deb-1) tampon=fgetc(ftblastx);
	paire[1]=tampon;
	paire[0]= ( (sens>0) ? X->AA(i,0) : X->AA(deb+fin-i-1,-1) );
	score= PROTMAT->VAL[PROTMAT->mot2indice(paire)];
	//	      printf("i:%d query-subject:  %c-%c score: %d\n",i,paire[0],paire[1],score);
      }
      GlobalBits=((REAL)bits)/((REAL)abs(fin-deb));
      TblastxNumber[phase][i]++;
      TblastxScore[phase][i]= ( (TblastxScore[phase][i]==0) ? score : Max ( (int)TblastxScore[phase][i],score)) ;
      //      TblastxScore[phase][i]= Max (TblastxScore[phase][i],GlobalBits);
    }
  }
  fclose(ftblastx);
  fprintf(stderr,"done\n");
  //	for(i=1495;i>1445;i-=3){
  //	for (i=1440;i<1500;i+=3) {
  //	  printf("i:%d  %c %c %c  %c %c %c\n",i,X->AA(i  ,1), X->AA(i+1,1), X->AA(i+2,1), X->AA(i  ,-1), X->AA(i+1,-1), X->AA(i+2,-1) );
  //	}
  
  //	  printf("i:%d  %c %c %c   %c\n",i,X->AA(i,-1),X->AA(i-1,-1),X->AA(i-2,-1),X->AA(i-3,-1)) ;
  
  if (PAR.getI("Output.graph")) Plot(X);
}

// --------------------------
//  GiveInfo Content TBlastX.
// --------------------------
void SensorHomology :: GiveInfo(DNASeq *X, int pos, DATA *d)
{
  int k;
  for (k = 0; k<6; k++) {
    if (TblastxNumber[k][pos]>0)
      d->ContentScore[k] += tblastxupdate(TblastxNumber[k][pos],TblastxScore[k][pos],TblastxP,TblastxB);
  }
}

// -----------------------------------------------
//  Convertit les phases 0-6 en 1 2 3 -1 -2 -3 0.
// -----------------------------------------------
int SensorHomology :: PhaseAdapt(char p)
{
  if (p >= 12) return 0;
  else if (p < 3) return (1+p);
  else if (p < 6) return (2-p);
  else if (p < 9) return (p-2);
  else return (5-p);
}

// -----------------------------------------------
//  Convertit les phases 1 2 3 -1 -2 -3 0 en 0-6.
// -----------------------------------------------
char SensorHomology :: ph06(char p)
{
  if (p == 0) return 6;
  else if (p > 0) return (p-1);
  else return 2-p;   
}
// ----------------------------
//  Plot Sensor information
// ----------------------------
void SensorHomology :: Plot(DNASeq *X)
{
  for (int pos =0; pos < X->SeqLen ; pos++){
    for(int j=0;j<6;j++) {
      if (TblastxNumber[j][pos]>0) {
	PlotBarI (pos,PhaseAdapt(j),0.6, TblastxNumber[j][pos], ((TblastxScore[j][pos]>0)? 8- (int) (Min(2,(int)TblastxScore[j][pos]/3)) : 3 ) );
      }
    }
  }
}

// --------------------------
//  Affichage des TBlastX.
// --------------------------
void SensorHomology :: PlotAt(int pos)
{
  //  if (PAR.getI("Output.graph")) {
  for(int j=0;j<6;j++) {
    if (TblastxNumber[j][pos]>0) {
      //	  PlotBarI (pos,PhaseAdapt(j),0.6, TblastxNumber[j][pos]/2 , 8- (int)(1.7*TblastxScore[j][pos]) );
      //  	    PlotBarI (pos,PhaseAdapt(j),0.6, TblastxNumber[j][pos], 6);//8- (int)(1.7*TblastxScore[j][pos]) );
      PlotBarI (pos,PhaseAdapt(j),0.6, TblastxNumber[j][pos], ((TblastxScore[j][pos]>0)? 8- (int) (Min(2,(int)TblastxScore[j][pos]/3)) : 3 ) );
    }
  }
}

// ------------------
//  Post analyse
// ------------------
void SensorHomology :: PostAnalyse(Prediction *pred)
{
}
