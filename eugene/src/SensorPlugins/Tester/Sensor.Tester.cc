/*****************************************************************************/
/*             Copyright (c) 2002 by INRA. All rights reserved.              */
/*                 Redistribution is not permitted without                   */
/*                 the express written permission of INRA.                   */
/*                     Mail : tschiex@toulouse.inra.fr                       */
/*---------------------------------------------------------------------------*/
/* File         : EuGeneTk/SensorPlugins/Tester/Sensor.Tester.cc             */
/* Description  : Sensor tester                                              */
/* Authors      : P.Bardou, S.Foissac, M.J.Cros, A.Moisan, T.Schiex          */
/* History      : May 2003         	   		                     */
/*****************************************************************************/

#include "Sensor.Tester.h"

extern Parameters PAR;

/*************************************************************
 **                      SensorTester                       **
 *************************************************************/
// ----------------------
//  Default constructor.
// ----------------------
SensorTester :: SensorTester (int n) : Sensor(n)
{
  char * sensorName = new char[FILENAME_MAX+1];
  char * paramKey   = new char[FILENAME_MAX+1];
  char * pluginsDir;
  char paramT[20] = "Tester.Test.";
  char paramI[20] = "Tester.Param.idx.";

  nbTest = 0;
  pluginsDir = PAR.getC("EuGene.PluginsDir");

  // Nombre de sensor à tester ?
  sprintf(paramKey,"%s%d",paramT,nbTest);
  while(PAR.count(paramKey)) {
    nbTest++;
    sprintf(paramKey,"%s%d",paramT,nbTest);
  }

  sensorL = new (SensorLoader *)[nbTest];
  sensor  = new (Sensor *)[nbTest];
  fp      = new (FILE *)[nbTest];
  source  = new (char *)[nbTest];
  
  for(int i=0; i<nbTest; i++) {
    sprintf(paramKey,"%s%d",paramT,i);
    source[i] = PAR.getC(paramKey);
    sprintf(sensorName,"%sSensor.%s.so",pluginsDir,source[i]);

    sensorL[i] =  new SensorLoader (sensorName);

    // On récupère le numéro d'instance
    sprintf(paramKey,"%s%d",paramI,i);
    int nbIdx = PAR.getI(paramKey);
    if (!sensorL[i]->LastError()) {
    fprintf(stderr," -Test.%d : Sensor.%.5s\t%d\n",i,source[i],nbIdx);
    sensor[i] = sensorL[i]->MakeSensor(nbIdx);
    }
    else {
      fprintf(stderr,"WARNING: ignored plugin (invalid or not found) : %s\n",
	      source[i]);
      exit(2);
    }
    
    char * outputFile = new char[FILENAME_MAX+1];
    sprintf(outputFile,"%stest.%s.gff",PAR.getC("Output.Prefix"),source[i]);
    // On verif outputFile n'existe pas
    if (!(fp[i] = fopen(outputFile, "r"))) {
      if (!(fp[i] = fopen(outputFile, "w"))) {
	fprintf(stderr, "cannot open %s output file\n", outputFile);
	exit(2);
      }
      else
	fprintf(fp[i],"SeqName\t Source\tFeature\t  Start\t    End\t  Score"
		"\t Strand\t  Frame\t    T/F\t  State\n");
    }
    else {
      fclose(fp[i]);
      fprintf(stderr, "WARNING: test output file \"%s\" exist\n", outputFile);
      exit(2);
    }
  }
}

// --------------------
//  Default destructor
// --------------------
SensorTester :: ~SensorTester ()
{
  for(int i=0; i<nbTest; i++) {
    delete sensorL[i];
    fclose(fp[i]);
  }
}

// -------------------
//  Init from new seq
// -------------------
void SensorTester :: Init (DNASeq *X)
{
  type = Type_Unknown;
  
  fprintf(stderr,"Test");
  for(int i=0; i<nbTest; i++)
    fprintf(stderr," %s",source[i]);
  fprintf(stderr," :\n");

  fprintf(stderr," Reading coordinates file......................");
  fflush(stderr);
  strcpy(seqName,PAR.getC("fstname"));
  strcat(seqName,".gff");
  ReadCoord(seqName);
  fprintf(stderr,"done\n");

  for(int i=0; i<nbTest; i++) {
    fprintf(stderr," ");
    sensor[i]->Init(X);
  }

  // Sequence name
  strcpy(seqName, BaseName(PAR.getC("fstname")));
  if (char * suffix = rindex(seqName,'.')) *suffix = 0;
}

// -----------------
//  Read coord file
// -----------------
void SensorTester :: ReadCoord(char name[FILENAME_MAX+1])
{
  // Lecture du fichier .gff (format gff) pour charger un objet
  // prediction "simplifié" (on distingue 6 états : IG, UTR, ExonF,
  // ExonR, IntronF et IntronR) -> but quel état réel à telle pos ?
  // Mots clefs pour le .gff : Features UTR5, UTR3, E.Init, E.Intr,
  // E.Term, E.Sngl

  FILE *fpCoord;
  int i;
  char *feature = new char[FILENAME_MAX];
  int  start,  end;
  char strand, frame;
  gene = new Prediction();

  if (!(fpCoord = fopen(name, "r"))) {
    fprintf(stderr, "cannot open gff file %s\n", name);
    exit(2);
  }
  int j=0;
  while (1) {
    j++;
    i = fscanf(fpCoord,"%*s %*s %s %d %d %*s %c %c",
	       feature, &start, &end, &strand, &frame);
    if (i < 5) {
      if (i==-1) {
	if(j==1)
	  fprintf(stderr,"WARNING: empty gff file !...");
      }
      else {
	fprintf(stderr, "\nError in gff file %s, line %d.\n", name, j);
	exit(2);
      }
    }
    else {
      if (j==1) {
	gene->add(start-1, InterGen5);
	if (strcmp(feature,"UTR5") == 0 || strcmp(feature,"UTR3") == 0)
	  gene->add(end, UTR5F);
	else if (strcmp(feature, "E.Init") == 0)
	  gene->add(end, ExonF1);
	else if (strcmp(feature, "E.Term") == 0)
	  gene->add(end, ExonR1);
	else if (strcmp(feature, "E.Sngl") == 0)
	  if (strand == '+') gene->add(end, ExonF1);
	  else               gene->add(end, ExonR1);
	else {
	  fprintf(stderr, "\n Error in gff file %s, line %d.\n", name, j);
	  fprintf(stderr, " WARNING :\n"
		  "   - Complete genes only in gff file.\n"
		  "   - Feature must be UTR5, UTR3, E.Init,"
		  " E.Intr, E.Term or E.Sngl.\n");
	  exit(2);
	}
      }
      else {
	if (strcmp(feature,"UTR5") == 0 || strcmp(feature,"UTR3") == 0)
	  gene->add(end, UTR5F);
	else if (strcmp(feature,"E.Init") == 0)
	  if (strand == '+') gene->add(end, ExonF1);
	  else {
	    gene->add(start-1, IntronR1);
	    gene->add(end, ExonR1);
	  }
	else if (strcmp(feature,"E.Term") == 0)
	  if (strand == '-') gene->add(end, ExonR1);
	  else {
	    gene->add(start-1, IntronF1);
	    gene->add(end, ExonF1);
	  }
	else if (strcmp(feature,"E.Sngl") == 0)
	  if (strand == '+') gene->add(end, ExonF1);
	  else               gene->add(end, ExonR1); 
	else if (strcmp(feature,"E.Intr") == 0) {
	  if (strand == '+') {
	    gene->add(start-1, IntronF1);
	    gene->add(end, ExonF1);
	  }
	  else {
	    gene->add(start-1, IntronR1);
	    gene->add(end, ExonR1);
	  }
	}
	else {
	  fprintf(stderr, "\n Error in gff file %s, line %d.\n", name, j);
	  fprintf(stderr, " %s : unknown feature (UTR5, UTR3, E.Init,"
		  " E.Intr, E.Term or E.Sngl).\n",feature);
	  exit(2);
	}
      }
    }
    char c = fgetc(fpCoord);
    while (c != '\n' && c != EOF) c=fgetc(fpCoord);
    if (fgetc(fpCoord) == EOF) break;
  }
  fclose(fpCoord);

  // Complete gene ?
  if (strcmp(feature, "E.Intr")  == 0  ||
      (strcmp(feature, "E.Init") == 0  &&  strand == '+')  ||
      (strcmp(feature, "E.Term") == 0  &&  strand == '-'))
    {
      fprintf(stderr, "\n Error in gff file %s, line %d.\n", name, j);
      fprintf(stderr, " WARNING : complete genes only in gff file.\n");
      exit(2);
    }
  
  // The prediction method have been develop for vector vPos & vESTMatch
  // describing the sequence form end to start...
  gene->reversePred();

  //printf("\n"); gene->print();
}

// ----------
//  GiveInfo
// ----------
void SensorTester :: GiveInfo (DNASeq *X, int pos, DATA *d)
{
  DATA Data;
  char * predSigType = new char[10];
  char * tf;
  char * truthState;

  for(int i=0; i<nbTest; i++) {
    for(int j=0; j<DATA::LastSigType; j++)
      Data.sig[j].Clear();
    sensor[i]->GiveInfo(X,pos,&Data);
    
    for(int j=DATA::Start; j<=DATA::Don; j++) {
      // Forward
      if (Data.sig[j].weight[Signal::Forward] != 0.0) {
	tf = SigType_TF(j, pos, &predSigType);

	truthState = State(pos);
	if (!strcmp(truthState, "ExonR") || !strcmp(truthState, "IntronR"))
	tf = "False";
	
	fprintf(fp[i],"%7.7s\t%7.7s\t%7s\t%7d\t      ."
		"\t%7.2f\t      +\t      .\t%7s\t%7s\n",
		seqName, source[i], predSigType, pos,
		Data.sig[j].weight[Signal::Forward],
		tf, truthState);
      }
      // Reverse
      if(Data.sig[j].weight[Signal::Reverse] != 0.0) {
	tf = SigType_TF(j, pos, &predSigType);

	truthState = State(pos);
	if (!strcmp(truthState, "ExonF") || !strcmp(truthState, "IntronF"))
	tf = "False";

	fprintf(fp[i],"%7.7s\t%7.7s\t%7s\t%7d\t      ."
		"\t%7.2f\t      -\t      .\t%7s\t%7s\n",
		seqName, source[i], predSigType, pos,
		Data.sig[j].weight[Signal::Reverse],
		tf, truthState);
      }
    }
  }
}

// -------------------------------------------
//  SigType_TF (Cf:enum SigType in SensorIF.h)
// -------------------------------------------
char* SensorTester :: SigType_TF(int i, int pos, char **sType)
{
  switch (i) {
  case 2:
    *sType = "Start";
    return gene->isStart(pos);
  case 3:
    *sType = "Stop";
    return gene->isStop(pos);
  case 4:
    *sType = "Acc";
    return gene->isAcc(pos);
  case 5:
    *sType = "Don";
    return gene->isDon(pos);
  default:
    *sType = "Unknown";
    return "- ? -";
  }
}

// -------
//  State
// -------
char* SensorTester :: State(int pos)
{
  switch (gene->getStateForPos(pos)) {
  case 0:
    return "ExonF";
  case 3:
    return "ExonR";
  case 6:
    return "IntronF";
  case 9:
    return "IntronR";
  case 13:
    return "UTR";
  case 12:
    return "IG";
  case -1:        // From the last element of "prediction"
    return "IG";  // to the end of the sequence...
  default:
    return ".";
  }
}

// -------------------------
//  Plot Sensor information
// -------------------------
void SensorTester :: Plot(DNASeq *X)
{
}

// --------------
//  Post analyse
// --------------
void SensorTester :: PostAnalyse(Prediction *pred)
{
}
