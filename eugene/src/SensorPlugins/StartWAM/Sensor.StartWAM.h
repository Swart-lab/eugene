#ifndef  SENSOR_StartWAM_INCLUDED
#define  SENSOR_StartWAM_INCLUDED

#include "../../EuGene/Sensor.h"
#include "../0_SensorTk/markov.h"

/*************************************************************
 **                 SensorStartWAM                    **
 *************************************************************/
class SensorStartWAM : public Sensor
{
 private:
  int order;    // order of the markov models in the Window Array Model
  double coef; // coefficient for the WAM score
  double pen; //  penality for the WAM score (score= coef * WAMscore - pen)
// merge flanking the consensus dinucleotide GT|AG (nbr of nt before and after)
  int beforestart,afterstart;
  int startsitelen;
  ChaineADN* ADN;
  TabChaine<ChaineADN,unsigned short int>* TPSTARTMOD; // True Positive Start site Model
  TabChaine<ChaineADN,unsigned short int>* FPSTARTMOD; // False         Start

 public:
  SensorStartWAM  (int);
  virtual ~SensorStartWAM   ();
  virtual void Init       (DNASeq *);
  virtual void GiveInfo   (DNASeq *, int, DATA *);
  virtual void Plot       (DNASeq *);
  virtual void PostAnalyse(Prediction *);
  void ReadModelFile (TabChaine<ChaineADN, unsigned short int> &model, char* filename);
  double startscoring (char* oligont,int order, int insitepos);
};

extern "C" SensorStartWAM * builder0(int n) {  return new SensorStartWAM(n); }

#endif
