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
    UTR3F = 11, UTR3R =12,
    IntronUTRF = 13, IntronUTRR = 14, 
    LastContentsType = 15};

  double  contents[LastContentsType];

  unsigned char EstMatch; 
};

// Type de sensor
enum TYPE_SENSOR {Type_Stop    = 1,    //               : Stop
		  Type_Start   = 2,    // Sensor SIGNAL : Start
		  Type_Acc     = 4,    //               : Acceptor
		  Type_Don     = 8,    //               : Donor
		  Type_TStop   = 16,   //               : Tr. Stop
		  Type_TStart  = 32,   //               : Tr. Start
		  Type_FS      = 64,   //               : FShift
		  Type_Content = 128,  //               : Contenu
		  Type_Any     = 255,  //               : Tous types
		  Type_None    = 0     //               : Type indefini
                  }; // Sensor de type inconnu (non encore initialisé)
#endif
