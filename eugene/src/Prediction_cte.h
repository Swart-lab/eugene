/** \file  Prediction_cte.h
 *  \brief  Les constantes utilisees par Prediction.h
 **/

#ifndef  __PREDICTION_CTE_H__
#define  __PREDICTION_CTE_H__


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
    IntronF2T = 30, IntronF3TG = 31, IntronF3TA = 32,
    IntronR3G = 33, IntronR3A  = 34, IntronR2AG = 35,
    InterGen  = 36,
    UTR5F     = 37, UTR3F      = 38,
    UTR5R     = 39, UTR3R      = 40,
    IntronU5F = 41, IntronU5R  = 42,
    IntronU3F = 43, IntronU3R  = 44,
    RnaF      = 45, RnaR       = 46,
    SnglF1F2  = 47, SnglF1F3   = 48, SnglF2F3  = 49,
    SnglR1R2  = 50, SnglR1R3   = 51, SnglR2R3  = 52,
    SnglF1R1  = 53, SnglF1R2   = 54, SnglF1R3  = 55,
    SnglF2R1  = 56, SnglF2R2   = 57, SnglF2R3  = 58, 
    SnglF3R1  = 59, SnglF3R2   = 60, SnglF3R3  = 61,
    UIRF      = 62, UIRR = 63,
    NbTracks  = 64
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

const int eukActiveTracksNb= 47;
const short int eukActiveTracks[eukActiveTracksNb]={
    InitF1,  InitF2,  InitF3 ,
    InitR1, InitR2,InitR3,
    SnglF1, SnglF2,SnglF3,
    SnglR1, SnglR2,SnglR3,
    IntrF1 ,IntrF2,IntrF3,
    IntrR1 ,IntrR2,IntrR3 ,
    TermF1,TermF2,TermF3,
    TermR1 ,TermR2,TermR3,
    IntronF1 , IntronF2 , IntronF3,
    IntronR1 , IntronR2 , IntronR3,
    IntronF2T , IntronF3TG, IntronF3TA,
    IntronR3G ,IntronR3A ,IntronR2AG,
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
    IntronR2AG, IntronR3A, IntronR3G,
    IntronF3TA, IntronF3TG, IntronF2T,
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

const short int ForwardTracks[(NbTracks-1)/2] =
{
    InitF1,InitF2,InitF3,
    SnglF1,SnglF2,SnglF3,
    IntrF1,IntrF2,IntrF3,
    TermF1,TermF2,TermF3,
    IntronF1,IntronF2,IntronF3,
    IntronF2T,IntronF3TG,IntronF3TA,
    UTR5F,UTR3F,
    IntronU5F,IntronU3F, RnaF,
    SnglF1F2, SnglF1F3, SnglF2F3,
    UIRF
};

const short int ReverseTracks[(NbTracks-1)/2] =
{
    InitR1,InitR2,InitR3,
    SnglR1,SnglR2,SnglR3,
    IntrR1,IntrR2,IntrR3,
    TermR1,TermR2,TermR3,
    IntronR1,IntronR2,IntronR3,
    IntronR3G,IntronR3A,IntronR2AG,
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
    SPLICED_TRANSCRIBED,SPLICED_TRANSCRIBED,SPLICED_TRANSCRIBED,
    SPLICED_TRANSCRIBED,SPLICED_TRANSCRIBED,SPLICED_TRANSCRIBED,
    SPLICED_TRANSCRIBED,SPLICED_TRANSCRIBED,SPLICED_TRANSCRIBED,
    SPLICED_TRANSCRIBED,SPLICED_TRANSCRIBED,SPLICED_TRANSCRIBED,
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
    FrameIntronF, FrameIntronF, FrameIntronF,
    FrameIntronR, FrameIntronR, FrameIntronR,
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
    Forward, Forward, Forward,
    Reverse, Reverse, Reverse,
    Forward, Forward, Forward,
    Reverse, Reverse, Reverse,   // Introns
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


#endif  //__PREDICTION_CTE_H__
