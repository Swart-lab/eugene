/** \file  Prediction_cte.h
 *  \brief  Les constantes utilisees par Prediction.h
 **/

#ifndef  __PREDICTION_CTE_H__
#define  __PREDICTION_CTE_H__

enum Tracks {
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
  NbTracks  = 45};



// Reverse the tracks
const enum Tracks ReverseIt[NbTracks] = {
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
  IntronU3R, IntronU3F
};

const short int UnorientedTracks[1] = {InterGen};

const short int ForwardTracks[(NbTracks-1)/2] = {
  InitF1,InitF2,InitF3,
  SnglF1,SnglF2,SnglF3,
  IntrF1,IntrF2,IntrF3,
  TermF1,TermF2,TermF3,
  IntronF1,IntronF2,IntronF3,
  IntronF2T,IntronF3TG,IntronF3TA,
  UTR5F,UTR3F, 
  IntronU5F,IntronU3F};

const short int ReverseTracks[(NbTracks-1)/2] = {
  InitR1,InitR2,InitR3,
  SnglR1,SnglR2,SnglR3,
  IntrR1,IntrR2,IntrR3,
  TermR1,TermR2,TermR3,
  IntronR1,IntronR2,IntronR3,
  IntronR3G,IntronR3A,IntronR2AG,
  UTR5R,UTR3R, 
  IntronU5R,IntronU3R};

// 0: untranscribed
// 1: transcribed but spliced
// 2: transcribed and in matured transcript but not translated
// 3: translated
const short int State2Status[NbTracks] = {3,3,3,3,3,3,
            3,3,3,3,3,3,
            3,3,3,3,3,3,
            3,3,3,3,3,3,
            1,1,1,1,1,1,
            1,1,1,1,1,1,
            0,
            2,2,2,2,
            1,1,1,1};


// 0 = interG, 1,2,3 = coding phase. 4 = intron. sign = strandness
const short int State2Phase[NbTracks] = {1,2,3,-1,-2,-3,
           1,2,3,-1,-2,-3,
           1,2,3,-1,-2,-3,
           1,2,3,-1,-2,-3,
            4, 4, 4,
           -4,-4,-4,
            4, 4, 4,
           -4,-4,-4,
           0,
           0, 0,0, 0,
           4,-4,4,-4};

// same with intron phase
const short int State2Frame[NbTracks] = {1,2,3,-1,-2,-3,
           1,2,3,-1,-2,-3,
           1,2,3,-1,-2,-3,
           1,2,3,-1,-2,-3,
           4,5,6,
           -4,-5,-6,
           5,6,6,
           -6,-6,-5,
           0,0,0,0,0,
           4,-4,4,-4};

#endif  //__PREDICTION_CTE_H__
