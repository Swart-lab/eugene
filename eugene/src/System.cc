#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#include "System.h"
#include <string.h>
#ifdef HAVE_STRINGS_H
#include <strings.h>
#endif


#include <time.h>
#include <string>
#include "Const.h"

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


// ------------------------------------------------------------------
// Provides the current date on the form jjMMMaaaa 
// where jj, aaaa are numbers and MMM caracters
// ------------------------------------------------------------------
void GetStrDate (char* d)
{

  char *m,*j,*a;
  m = new char[4]; j = new char[3]; a = new char[5];
  time_t t=time(0);

  strcpy(d,ctime(&t));
  sscanf(d, "%*s %s %s %*s %s", m,j,a);
  strcpy(d,j);strcat(d,m);strcat(d,a);

  delete [] m; delete [] j; delete [] a;
}
