typedef struct UTIL *ptUTIL ;

enum CHECK {check,nocheck};

typedef struct UTIL
{
  int tab;
  int   n1;
  int   n2;
  double  delta;
  char* rais;
  enum CHECK check;
  ptUTIL suiv;
} UTIL;

extern int Utilisateur(char *nom_fich,ptUTIL *a, ptUTIL *b);
void Util(int i, ptUTIL ut, DATA *d);
void WriteUtils(ptUTIL,FILE *);
