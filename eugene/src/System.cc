#include "System.h"
#include <string.h>

// ------------------------------------------------------------------
// BASENAME: returns a pointer to the filename, w/o any
// leading prefix
// ------------------------------------------------------------------
char * BaseName(char *path)
{
  char *lead = rindex(path,'/');
  
  return ((lead != NULL)  ? lead+1 : path);
}
// ------------------------------------------------------------------
// Function to open a file by using environment variables.  We first
// search in local directory and if not found in directory defdir
// ------------------------------------------------------------------

FILE *FileOpen(const char *defdir, const char *filename, const char *mode)
{
  FILE  *fp;
  char buffer[FILENAME_MAX];
  
  if ((fp = fopen(filename, mode)))
    return fp;
  
  if (defdir) {
    strcat(strcat(strcpy(buffer, defdir), "/"), filename);
    fp = fopen (buffer, mode);
  }
  
  if  (!fp) {
    fprintf (stderr, "ERROR:  Could not open file %s \n", filename);
    exit (1);
  }  
  return  fp;
}
// ------------------------------------------------------------------
// Malloc sur
// ------------------------------------------------------------------

void *  Safe_malloc  (size_t Len)     
{
  void  * P;
  
  
  if  ((P = malloc (Len)) == NULL)
    {
      perror("malloc");
      exit (1);
    }
  return  P;
}

// ------------------------------------------------------------------
// Realloc sur
// ------------------------------------------------------------------

void *  Safe_realloc  (void * Q, size_t Len)
{
  void  * P;
  
  
  if  ((P = realloc (Q, Len)) == NULL)
    {
      perror("realloc");
      exit (1);
    }  
  return  P;
}

// -----------------------------------------------
// Gestion des arguments (integer)
// -----------------------------------------------

int GetIArg(char *arg, int *value, int dft)
{
    int varg;
    int ok;

    ok = (sscanf(arg, "%d", &varg) == 1) && (varg >= 0);

    *value = (ok ? varg : dft);
 
    return ok;
}
// -----------------------------------------------
// Gestion des arguments (double)
// -----------------------------------------------

int GetDArg(char *arg, double *value, double dft)
{
    double varg;
    int ok;

    ok = (sscanf(arg, "%lf", &varg) == 1) && (varg >= 0);

    *value = (ok ? -varg : dft);
 
    return ok;
}

