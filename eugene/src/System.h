#ifndef SYSTEM_H
#define SYSTEM_H

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#define myassert(x) assert(x)
//#define myassert(x)

#if defined (__sun)
int isinf(double x);
#endif

char * BaseName(char *path);
int ProbeFile(char *path);


template<class T>
inline T Min(T x, T y)
{
  if (x < y) return x;
  else
    return y;
}

template<class T>
inline T Max(T x, T y)
{
  if (x < y) return y;
  else
    return x;
}

FILE *FileOpen(const char *defdir, const char *filename, const char *mode);
void *  Safe_malloc  (size_t Len);
void *  Safe_realloc  (void * Q, size_t Len);
double cpuTime();
void GetStrDate (char* d);

#endif
