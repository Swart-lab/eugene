#ifndef  SENSOR_API_H_INCLUDED
#define  SENSOR_API_H_INCLUDED

class Signal
{
 public:

  enum Edge {Forward = 0, Reverse = 1, ForwardNo = 2, ReverseNo = 3, LastEdge = 4};
  double weight[LastEdge];

  inline void Clear() { for (int i = 0; i < LastEdge; i++) weight[i] = 0.0;}
  inline bool IsSet(int i) { return ((weight[i] != 0.0) || (weight[i+2] != 0.0)); }
  inline void SetToDefault() { 
    if (!IsSet(Forward)) weight[Forward] = NINFINITY;
    if (!IsSet(Reverse)) weight[Reverse] = NINFINITY;
  }
};

class DATA
{
 public:
  // Signal scores
  enum SigType {tStart = 0, tStop = 1, Start = 2, Stop = 3,
		Acc = 4, Don = 5, Ins = 6, Del =7, LastSigType = 8};
  Signal sig[LastSigType];

  // Contents scores 
  enum ContentsType {
    ExonF1 =0, ExonF2 =1, ExonF3 =2,
    ExonR1 = 3, ExonR2 =4, ExonR3 =5,
    IntronF = 6, IntronR = 7,
    InterG = 8,
    UTR5F = 9, UTR5R = 10,
    UTR3F = 11, UTR3R =12, LastContentsType = 13};

  double  contents[LastContentsType];

  // WARNING temporaire : EST -> on est dans intron
  unsigned char ESTMATCH_TMP; 
};

// Type de sensor
enum TYPE_SENSOR {Type_Stop,     //               : Stop
		  Type_Start,    // Sensor SIGNAL : Start
		  Type_Splice,   //               : Splice
		  Type_Content,  // Sensor CONTENT
		  Type_Multiple, // Sensor porteur de multiple info (Ex : User)
		  Type_Unknown}; // Sensor de type inconnu (non encore initialisé)


#endif
