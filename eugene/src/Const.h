#ifndef  CONST_H_INCLUDED
#define  CONST_H_INCLUDED
#include <math.h>
#include <stdlib.h>

#define EugDir getenv("EUGENEDIR")

//#define REAL double
#define REAL float

const int FASTA_Len = 50;
const int MAX_LINE = 300;
const int TRUE  = 1;
const int FALSE = 0;

const double NINFINITY   = log(0.0);
const REAL DontCrossStop = NINFINITY;

typedef struct DATA {
  REAL Stop[2];
  REAL Start[2];
  REAL Acc[2];
  REAL Don[2];
  REAL ContentScore[13];
  unsigned char ESTMATCH_TMP; // WARNING temporaire : EST -> on est dans intron
}DATA;

// Type de sensor
enum TYPE_SENSOR {Type_Stop,     //               : Stop
		  Type_Start,    // Sensor SIGNAL : Start
		  Type_Splice,   //               : Splice
		  Type_Content,  // Sensor CONTENT
		  Type_Multiple, // Sensor porteur de multiple info (Ex : User)
		  Type_Unknown}; // Sensor de type inconnu (non encore initialisé)

// Les etats
const int ExonF1 = 0;
const int ExonF2 = 1;
const int ExonF3 = 2;
const int ExonR1 = 3;
const int ExonR2 = 4;
const int ExonR3 = 5;

const int IntronF1 = 6;
const int IntronF2 = 7;
const int IntronF3 = 8;
const int IntronR1 = 9;
const int IntronR2 = 10;
const int IntronR3 = 11;

const int InterGen5 = 12;
const int InterGen3 = 17;

const int UTR5F = 13;
const int UTR5R = 15;
const int UTR3F = 14;
const int UTR3R = 16;

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
