#ifndef SYSTEM_H
#define SYSTEM_H

#define  TRUE  1
#define  FALSE  0

#define MIN(x, y) (((x) < (y)) ? (x) : (y))
#define MAX(x, y) (((x) < (y)) ? (y) : (x))


FILE *  File_Open  (const char * Filename, const char * Mode)

{
  FILE  *  fp;
  
  if  ((fp = fopen (Filename, Mode)) == NULL)
    {
      fprintf (stderr, "ERROR:  Could not open file  %s \n", Filename);
      exit (1);
    }  
  return  fp;
  }

// Malloc sur

void *  Safe_malloc  (size_t Len)     
{
  void  * P;
  
  
  if  ((P = malloc (Len)) == NULL)
    {
      fprintf (stderr, "ERROR:  malloc failed\n");
      exit (1);
    }
  return  P;
}

// Realloc sur

void *  Safe_realloc  (void * Q, size_t Len)
{
  void  * P;
  
  
  if  ((P = realloc (Q, Len)) == NULL)
    {
      fprintf (stderr, "ERROR:  realloc failed\n");
      exit (1);
    }  
  return  P;
}

#endif
