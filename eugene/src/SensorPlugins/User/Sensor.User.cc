#include "SensorUser.h"

/*************************************************************
 **                        SensorUser                       **
 *************************************************************/

extern Parameters PAR;

// ----------------------
//  Default constructor.
// ----------------------
SensorUser :: SensorUser ()
{
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

  type = Type_Multiple;

  fprintf(stderr,"Loading user data...");
  
  strcpy(tempname, PAR.getC("fstname"));
  strcat(tempname, ".user");
  errflag = Utilisateur(tempname);   //prise en compte de donnees utilisateur
    
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
  UserInfoList = SignalUser;
  Util(pos, UserInfoList, d);
  
  UserInfoList = ContentsUser;
  Util(pos, UserInfoList, d);
}
