#include "Sensor.User.h"
#include "yacc.tab.h"

/*************************************************************
 **                        SensorUser                       **
 *************************************************************/

extern Parameters PAR;

// ----------------------
//  Default constructor.
// ----------------------
SensorUser :: SensorUser (int n, DNASeq *X) : Sensor(n)
{
  type = Type_Multiple;
}

// ----------------------
//  Default destructor.
// ----------------------
SensorUser :: ~SensorUser ()
{
}

// ----------------------
//  Init user.
// ----------------------
void SensorUser :: Init (DNASeq *X)
{
  char tempname[FILENAME_MAX+1];
  int errflag;

  fprintf(stderr,"Reading user data............");
  
  strcpy(tempname, PAR.getC("fstname"));
  strcat(tempname, ".user");
  errflag = Utilisateur(tempname, &Signals, &Contents);   //prise en compte de donnees utilisateur
    
  if (errflag) {
    fprintf(stderr,"none found\n");
  }
  else fprintf(stderr,"done\n");
}

// -----------------------
//  GiveInfo signal user.
// -----------------------
void SensorUser :: GiveInfo (DNASeq *X, int pos, DATA *d)
{
  Util(pos, Signals, d);
  Util(pos, Contents, d);
}

// ----------------------------
//  Plot Sensor information
// ----------------------------
void SensorUser :: Plot(DNASeq *X)
{
}

// ------------------
//  Post analyse
// ------------------
void SensorUser :: PostAnalyse(Prediction *pred)
{
}
