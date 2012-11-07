/** \file  Prediction_cte.h
 *  \brief  Les constantes utilisees par Prediction.h
 **/

#ifndef  __PREDICTION_CTE_H__
#define  __PREDICTION_CTE_H__

#define DEV

enum Tracks
{
    InitF1    = 0,  InitF2     = 1,  InitF3     = 2,
    InitR1    = 3,  InitR2     = 4,  InitR3     = 5,
    SnglF1    = 6,  SnglF2     = 7,  SnglF3     = 8,
    SnglR1    = 9,  SnglR2     = 10, SnglR3     = 11,
    IntrF1    = 12, IntrF2     = 13, IntrF3     = 14,
    IntrR1    = 15, IntrR2     = 16, IntrR3     = 17,
    TermF1    = 18, TermF2     = 19, TermF3     = 20,
    TermR1    = 21, TermR2     = 22, TermR3     = 23,
    IntronF1  = 24, IntronF2   = 25, IntronF3   = 26,
    IntronR1  = 27, IntronR2   = 28, IntronR3   = 29,
    IntronF2T = 30, IntronF2A  = 31, IntronF3TG = 32,                   // Txx,   Axx,   TG[A]
    IntronF3TA = 33, IntronF3TC = 34, IntronF3AG = 35,                  // TA[GA],TC[A], AG[AG]
    IntronR3G  = 36, IntronR3A  = 37, IntronR2GA = 38,                  // xxG,   xxA,   [TA]GA
    IntronR2AG = 39, IntronR2AA = 40, IntronR2CA = 41, IntronR2GG = 42, // [T]AG, [T]AA, [T]CA, [A]GG
    InterGen  = 43,
    UTR5F     = 44, UTR3F      = 45,
    UTR5R     = 46, UTR3R      = 47,
    IntronU5F = 48, IntronU5R  = 49,
    IntronU3F = 50, IntronU3R  = 51,
    RnaF      = 52, RnaR       = 53,
    SnglF1F2  = 54, SnglF1F3   = 55, SnglF2F3  = 56,
    SnglR1R2  = 57, SnglR1R3   = 58, SnglR2R3  = 59,
    SnglF1R1  = 60, SnglF1R2   = 61, SnglF1R3  = 62,
    SnglF2R1  = 63, SnglF2R2   = 64, SnglF2R3  = 65, 
    SnglF3R1  = 66, SnglF3R2   = 67, SnglF3R3  = 68,
    UIRF      = 69, UIRR = 70,
    NbTracks  =  71 
};

const int proActiveTracksNb= 30;
const short int proActiveTracks[proActiveTracksNb]={
SnglF1, SnglF2, SnglF3,
SnglR1, SnglR2, SnglR3,
InterGen,
UTR5F, UTR3F,
UTR5R, UTR3R,
RnaF, RnaR,
SnglF1F2, SnglF1F3, SnglF2F3,
SnglR1R2, SnglR1R3, SnglR2R3,
SnglF1R1, SnglF1R2, SnglF1R3,
SnglF2R1, SnglF2R2, SnglF2R3,
SnglF3R1, SnglF3R2, SnglF3R3,
UIRF, UIRR
};

const int proForwardActiveTracksNb= 11;
const short int proForwardActiveTracks[proForwardActiveTracksNb]={
SnglF1, SnglF2, SnglF3,
InterGen,
UTR5F, UTR3F,
RnaF,
SnglF1F2, SnglF1F3, SnglF2F3,
UIRF
};

const int proReverseActiveTracksNb= 11;
const short int proReverseActiveTracks[proReverseActiveTracksNb]={
SnglR1, SnglR2, SnglR3,
InterGen,
UTR5R, UTR3R,
RnaR,
SnglR1R2, SnglR1R3, SnglR2R3,
UIRR
};

const int eukActiveTracksNb= 54;
const short int eukActiveTracks[eukActiveTracksNb]={
    InitF1,  InitF2,  InitF3 ,
    InitR1, InitR2,InitR3,
    SnglF1, SnglF2,SnglF3,
    SnglR1, SnglR2,SnglR3,
    IntrF1 ,IntrF2,IntrF3,
    IntrR1 ,IntrR2,IntrR3 ,
    TermF1,TermF2,TermF3,
    TermR1 ,TermR2,TermR3,
    IntronF1, IntronF2 , IntronF3,
    IntronR1, IntronR2, IntronR3,
    IntronF2T,  IntronF2A,  IntronF3TG,                   // Txx,   Axx,   TG[A]
    IntronF3TA, IntronF3TC, IntronF3AG,                 // TA[GA],TC[A], AG[AG]
    IntronR3G,  IntronR3A,  IntronR2GA,                  // xxG,   xxA,   [TA]GA
    IntronR2AG, IntronR2AA, IntronR2CA, IntronR2GG,
    InterGen,
    UTR5F, UTR3F,
    UTR5R , UTR3R ,
    IntronU5F,IntronU5R,
    IntronU3F,IntronU3R,
    RnaF, RnaR
};



// Reverse the tracks
const enum Tracks ReverseIt[NbTracks] =
{
    InitR1, InitR2, InitR3,
    InitF1, InitF2, InitF3,
    SnglR1, SnglR2, SnglR3,
    SnglF1, SnglF2, SnglF3,
    IntrR1, IntrR2, IntrR3,
    IntrF1, IntrF2, IntrF3,
    TermR1, TermR2, TermR3,
    TermF1, TermF2, TermF3,
    IntronR1, IntronR2, IntronR3,
    IntronF1,IntronF2, IntronF3,
    IntronR2AG, IntronR2GA,  IntronR3A,                   // Txx,   Axx,   TG[A]
    IntronR3G,  IntronR3A,   IntronR3G,                 // TA[GA],TC[A], AG[AG]
    IntronF3TA, IntronF3TA,  IntronF2T,                  // xxG,   xxA,   [TA]GA
    IntronF2T,  IntronF2T,   IntronF2T, IntronF2A,
    InterGen,
    UTR5R,UTR3R,
    UTR5F, UTR3F,
    IntronU5R, IntronU5F,
    IntronU3R, IntronU3F,
    RnaR,      RnaF,
    SnglR1R2, SnglR1R3, SnglR2R3,
    SnglF1F2, SnglF1F3, SnglF2F3,
    SnglF1R1, SnglF2R1, SnglF3R1,
    SnglF1R2, SnglF2R2, SnglF3R2,
    SnglF1R3, SnglF2R3, SnglF3R3,
    UIRR, UIRF
};

const short int UnorientedTracks[1] = {InterGen};

const short int ForwardTracks[30] =
{
    InitF1,InitF2,InitF3,
    SnglF1,SnglF2,SnglF3,
    IntrF1,IntrF2,IntrF3,
    TermF1,TermF2,TermF3,
    IntronF1,IntronF2,IntronF3,
    IntronF2T,  IntronF2A,  IntronF3TG,
    IntronF3TA, IntronF3TC, IntronF3AG,  
    UTR5F,UTR3F,
    IntronU5F,IntronU3F, RnaF,
    SnglF1F2, SnglF1F3, SnglF2F3,
    UIRF
};

const short int ReverseTracks[31] =
{
    InitR1,InitR2,InitR3,
    SnglR1,SnglR2,SnglR3,
    IntrR1,IntrR2,IntrR3,
    TermR1,TermR2,TermR3,
    IntronR1,IntronR2,IntronR3,
    IntronR3G,IntronR3A,IntronR2GA,             
    IntronR2AG,IntronR2AA,IntronR2CA,IntronR2GG,
    UTR5R,UTR3R,
    IntronU5R,IntronU3R, RnaR,
    SnglR1R2, SnglR1R3, SnglR2R3,
    UIRR
};

// 0: untranscribed
// 1: transcribed but spliced
// 2: transcribed and in matured transcript but not translated
// 3: translated
enum Status
{
    UNTRANSCRIBED       = 0,
    SPLICED_TRANSCRIBED = 1,
    UNTRANSLATED        = 2,
    TRANSLATED          = 3,
    RNA                 = 4
};


const short int State2Status[NbTracks] =
{
    TRANSLATED,TRANSLATED,TRANSLATED,
    TRANSLATED,TRANSLATED,TRANSLATED,
    TRANSLATED,TRANSLATED,TRANSLATED,
    TRANSLATED,TRANSLATED,TRANSLATED,
    TRANSLATED,TRANSLATED,TRANSLATED,
    TRANSLATED,TRANSLATED,TRANSLATED,
    TRANSLATED,TRANSLATED,TRANSLATED,
    TRANSLATED,TRANSLATED,TRANSLATED,
    SPLICED_TRANSCRIBED,SPLICED_TRANSCRIBED,SPLICED_TRANSCRIBED, //introns
    SPLICED_TRANSCRIBED,SPLICED_TRANSCRIBED,SPLICED_TRANSCRIBED,
    SPLICED_TRANSCRIBED,SPLICED_TRANSCRIBED,SPLICED_TRANSCRIBED, // special introns
    SPLICED_TRANSCRIBED,SPLICED_TRANSCRIBED,SPLICED_TRANSCRIBED,
    SPLICED_TRANSCRIBED,SPLICED_TRANSCRIBED,SPLICED_TRANSCRIBED,
    SPLICED_TRANSCRIBED,SPLICED_TRANSCRIBED,SPLICED_TRANSCRIBED,SPLICED_TRANSCRIBED,
    UNTRANSCRIBED,
    UNTRANSLATED,UNTRANSLATED,
    UNTRANSLATED,UNTRANSLATED,
    SPLICED_TRANSCRIBED,SPLICED_TRANSCRIBED,
    SPLICED_TRANSCRIBED,SPLICED_TRANSCRIBED,
    RNA, RNA,
    TRANSLATED, TRANSLATED, TRANSLATED,
    TRANSLATED, TRANSLATED, TRANSLATED,
    TRANSLATED, TRANSLATED, TRANSLATED,
    TRANSLATED, TRANSLATED, TRANSLATED,
    TRANSLATED, TRANSLATED, TRANSLATED,
    UNTRANSLATED, UNTRANSLATED
};


enum Frame
{
    Frame1F      =  1, Frame2F      =  2, Frame3F =  3,
    Frame1R      = -1, Frame2R      = -2, Frame3R = -3,
    FrameIntronF =  4, FrameIntronR = -4, FrameIG =  0
};



// 0 = interG, 1,2,3 = coding phase. 4 = intron or rna. sign = strandness
// This array is just necessary for coding exons and introns
// ncRNA frame is FrameIG, to allow to display rna to the Intergenic track in graphic output
const short int State2Frame[NbTracks] =
{
    Frame1F, Frame2F, Frame3F,
    Frame1R, Frame2R, Frame3R,
    Frame1F, Frame2F, Frame3F,
    Frame1R, Frame2R, Frame3R,
    Frame1F, Frame2F, Frame3F,
    Frame1R, Frame2R, Frame3R,
    Frame1F, Frame2F, Frame3F,
    Frame1R, Frame2R, Frame3R,
    FrameIntronF, FrameIntronF, FrameIntronF,
    FrameIntronR, FrameIntronR, FrameIntronR,
    FrameIntronF, FrameIntronF, FrameIntronF, //special introns
    FrameIntronF, FrameIntronF, FrameIntronF,
    FrameIntronR, FrameIntronR, FrameIntronR,
    FrameIntronR, FrameIntronR, FrameIntronR, FrameIntronR,
    FrameIG,
    FrameIG, FrameIG,
    FrameIG, FrameIG,
    FrameIntronF, FrameIntronR,
    FrameIntronF, FrameIntronR,
    FrameIG, FrameIG,			// rna tracks.
    FrameIG, FrameIG, FrameIG,  // never used so initialized with FrameIG value
    FrameIG, FrameIG, FrameIG,  // never used so initialized with FrameIG value
    FrameIG, FrameIG, FrameIG,  // never used so initialized with FrameIG value
    FrameIG, FrameIG, FrameIG,  // never used so initialized with FrameIG value
    FrameIG, FrameIG, FrameIG,  // never used so initialized with FrameIG value
    FrameIG, FrameIG,
};


enum enumStrand
{
	Forward = 0, Reverse = 1, Forward_Reverse = 2, No_Strand = 3
};

const short int State2Strand[NbTracks] =
{
    Forward, Forward, Forward,
    Reverse, Reverse, Reverse,
    Forward, Forward, Forward,
    Reverse, Reverse, Reverse,
    Forward, Forward, Forward,
    Reverse, Reverse, Reverse,
    Forward, Forward, Forward,
    Reverse, Reverse, Reverse,   // Exons
    Forward, Forward, Forward,  // introns
    Reverse, Reverse, Reverse,
    Forward, Forward, Forward,   // special introns
    Forward, Forward, Forward,
    Reverse, Reverse, Reverse,   
    Reverse, Reverse, Reverse, Reverse,
    No_Strand,                   // IG
    Forward, Forward,
    Reverse, Reverse,            // UTR
    Forward, Reverse,
    Forward, Reverse,            // Intron
    Forward, Reverse,            // Rna
    Forward, Forward, Forward,   // Bicoding
    Reverse, Reverse, Reverse,   // Bicoding
    Forward_Reverse, Forward_Reverse, Forward_Reverse, // Bicoding
    Forward_Reverse, Forward_Reverse, Forward_Reverse, // Bicoding
    Forward_Reverse, Forward_Reverse, Forward_Reverse, // Bicoding
    Forward, Reverse,
};

const int AllowedStopNb = 6;
const char* const AllowedStop[AllowedStopNb] = {"tag", "tga", "taa", "aga", "agg", "tca"};

#endif  //__PREDICTION_CTE_H__
