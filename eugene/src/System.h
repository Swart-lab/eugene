#ifndef SYSTEM_H
#define SYSTEM_H

#define  TRUE  1
#define  FALSE 0

#define MIN(x, y) (((x) < (y)) ? (x) : (y))
#define MAX(x, y) (((x) < (y)) ? (y) : (x))

inline FILE *FileOpen(const char *defdir, const char *filename, const char *mode)
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

// Malloc sur
inline void *  MyMalloc(size_t Len)     
{
  void  *Ptr;
  
  if  (!(Ptr = malloc (Len))) {
    fprintf (stderr, "ERROR:  malloc failed\n");
    exit (1);
  }
  return  Ptr;
}

// Realloc sur

inline void *  MyRealloc(void *Or, size_t Len)
{
  void  *Ptr;
  
  if  (!(Ptr = realloc (Or, Len))) {
    fprintf (stderr, "ERROR:  realloc failed\n");
    exit (1);
  }  
  return Ptr;
}


#endif
