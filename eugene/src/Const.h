#ifndef  CONST_H_INCLUDED
#define  CONST_H_INCLUDED
#include <math.h>
#include <stdlib.h>

#define EugDir getenv("EUGENEDIR")

const int FASTA_Len = 50;
const int MAX_LINE = 300;
const int TRUE  = 1;
const int FALSE = 0;
const double NINFINITY   = log(0.0);

// Les Hits EST
const unsigned char HitForward    = 0x1;
const unsigned char MLeftForward  = 0x2;
const unsigned char GapForward    = 0x4;
const unsigned char MRightForward = 0x8; 

// Shift to go from Hits to ...
const  unsigned int HitToMLeft  = 1;
const  unsigned int HitToGap    = 2;
const  unsigned int HitToMRight = 3;

const unsigned char HitReverse    = 0x10;
const unsigned char MLeftReverse  = 0x20;
const unsigned char GapReverse    = 0x40;
const unsigned char MRightReverse = 0x80; 

const unsigned char Hit           = HitForward    | HitReverse;
const unsigned char MLeft         = MLeftForward  | MLeftReverse;
const unsigned char MForward      = MLeftForward  | MRightForward;
const unsigned char MReverse      = MLeftReverse  | MRightReverse;
const unsigned char Gap           = GapForward    | GapReverse;
const unsigned char MRight        = MRightForward | MRightReverse;
const unsigned char Margin        = MRight        | MLeft;

const unsigned char NotAHit       = Margin | Gap;

#endif
