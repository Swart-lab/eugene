#ifndef  PENALTYD_H_INCLUDED
#define  PENALTYD_H_INCLUDED

// A class to represent length penalties. The penalty can be arbitrary
// up to a given length and then linearly interpolated for further
// length.  Considering prob. distributions this means we have an
// exponential tail. INFINITE costs are allowed.

#include <vector>
#include <math.h>
#include <cstdio>
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

class PenaltyDist
{
  friend class Track;
 private:

  int MinLen;
  int MaxLen;
  std::vector<double> Distribution;
  std::vector<double> MinDist;
  std::vector<double> Delta;
  double FinalSlope;
  double MaxDelta;

 public:
  PenaltyDist();
  ~PenaltyDist();

  void   UpdateMinDist();
  void   UpdateMaxDelta ();

  inline double PenaltyDist :: GetDelta(int len) {
    if (len < (int)Delta.size()) return Delta[len];
    return MaxDelta;
  }
  double MinPen(int len);
  void   LoadPenaltyDist(char * FileName);
  void   PlotExp(FILE* h,int a,int b);
  void   Dump(FILE* h);
  double Integrate();
  void   Normalize();
  inline double operator [] (int len)
    {
      if (len < MinLen) return -log(0.0);
      if (len - MinLen < (int)Distribution.size())
	return Distribution[len-MinLen];
      return Distribution.back()+(len-MinLen-Distribution.size()+1)*FinalSlope;
    };
};


#endif
