/*****************************************************************************/
/*             Copyright (c) 2002 by INRA. All rights reserved.              */
/*                 Redistribution is not permitted without                   */
/*                 the express written permission of INRA.                   */
/*                     Mail : tschiex@toulouse.inra.fr                       */
/*---------------------------------------------------------------------------*/
/* File         : EuGeneTk/SensorPlugins/Plotter/Sensor.Plotter.cc           */
/* Description  : Sensor Plotter                                             */
/* Authors      : P.Bardou, S.Foissac, M.J.Cros, A.Moisan, T.Schiex          */
/* History      : May 2003                                                   */
/*****************************************************************************/

#include "Sensor.Plotter.h"

/*************************************************************
 **                       SensorPlotter
 *************************************************************/

extern Parameters PAR;

// ----------------------
// Default constructor.
// ----------------------
SensorPlotter :: SensorPlotter (int n, DNASeq *X) : Sensor(n)
{
  type = Type_Multiple;

  // Save parameters to limit the map access number
  window  = PAR.getI("Output.window");
  plotGC  = PAR.getI("Plotter.GC");
  plotGC3 = PAR.getI("Plotter.GC3");
  plotAT  = PAR.getI("Plotter.A|T/A+T");
}

// ----------------------
//  Default destructor.
// ----------------------
SensorPlotter :: ~SensorPlotter ()
{
}

// ----------------------
//  Init.
// ----------------------
void SensorPlotter :: Init (DNASeq *X)
{
  if (PAR.getI("Output.graph")) Plot(X);
}

// ------------------------
//  GiveInfo.
// ------------------------
void SensorPlotter :: GiveInfo (DNASeq *X, int pos, DATA *d)
{
}

// ----------------------------
//  Plot Sensor information.
// ----------------------------
void SensorPlotter :: Plot(DNASeq *X)
{
  double GC[3], pGC[3], A_AT, pA_AT, T_AT, pT_AT;
  unsigned int Width, pWidth;
  unsigned int GC1, AT1, A, T;
  int i;
  
  GC[0] = GC[1] = GC[2] = 0;
  A_AT = 0;
  T_AT = 0;
  
  for (i = 0; i < window/2; i++) {
    GC1 = (((*X)(i,0) & CodeG) != 0) + (((*X)(i,0) & CodeC) != 0);
    AT1 = (((*X)(i,0) & CodeA) != 0) + (((*X)(i,0) & CodeT) != 0);
    GC[i%3] += ((double)GC1/(GC1+AT1));
    
    A = (((*X)(i,0) & CodeA) != 0);
    T = (((*X)(i,0) & CodeT) != 0);
    if (AT1 != 0) {
      A_AT += ((double)A/AT1);
      T_AT += ((double)T/AT1);
    }
  }
  
  Width = window/2;
  
  for (i = 0; i < X->SeqLen; i++) {
    pGC[0] = GC[0];
    pGC[1] = GC[1];
    pGC[2] = GC[2];
    pA_AT  = A_AT;
    pT_AT  = T_AT;

    pWidth = Width;
    
    if (i-window/2 >= 0) {
      GC1 = (((*X)(i-window/2,0) & CodeG) != 0) +
	(((*X)(i-window/2,0) & CodeC) != 0);
      AT1 = (((*X)(i-window/2,0) & CodeA) != 0) +
	(((*X)(i-window/2,0) & CodeT) != 0);
      GC[(i-window/2)%3] -= ((double)GC1/(GC1+AT1));
      Width--;

      A = (((*X)(i-window/2,0) & CodeA) != 0);
      T = (((*X)(i-window/2,0) & CodeT) != 0);
      if (AT1 != 0) {
	A_AT -= ((double)A/AT1);
	T_AT -= ((double)T/AT1);
      }
    }

    if (i+window/2 < X->SeqLen) {
      GC1 = (((*X)(i+window/2,0) & CodeG) != 0) +
	(((*X)(i+window/2,0) & CodeC) != 0);
      AT1 = (((*X)(i+window/2,0) & CodeA) != 0) +
	(((*X)(i+window/2,0) & CodeT) != 0);
      GC[(i+window/2)%3] += ((double)GC1/(GC1+AT1));
      Width++;
      
      A = (((*X)(i+window/2,0) & CodeA) != 0);      
      T = (((*X)(i+window/2,0) & CodeT) != 0);
      if (AT1 != 0) {
	A_AT += ((double)A/AT1);
	T_AT += ((double)T/AT1);
      }
    }

    if (plotGC)
      PlotLine(((i == 0) ? 0 : i-1),i, 0, 0,
               ((double)(pGC[0]+pGC[1]+pGC[2])/ pWidth)*1.5 - 0.25,
               ((double)(GC[0] +GC[1] +GC[2]) / Width) *1.5 - 0.25, 5);
    
    if (plotGC3) {
      PlotLine(((i == 0) ? 0 : i-1),i, 1, 1,
               (((double)pGC[2]*3) / pWidth),
               (((double)GC[2] *3) / Width), 5);
      
      PlotLine(((i == 0) ? 0 : i-1),i, 2, 2,
               (((double)pGC[0]*3) / pWidth),
               (((double)GC[0] *3) / Width), 5);
      
      PlotLine(((i == 0) ? 0 : i-1),i, 3, 3,
               (((double)pGC[1]*3) / pWidth),
               (((double)GC[1] *3) / Width), 5);
      
      PlotLine(((i == 0) ? 0 : i-1),i, -1, -1,
               (((double)pGC[X->SeqLen%3]*3) / pWidth),
               (((double)GC[X->SeqLen%3] *3) / Width), 5);
      
      PlotLine(((i == 0) ? 0 : i-1),i, -2, -2,
               (((double)pGC[(X->SeqLen-1)%3]*3) / pWidth),
               (((double)GC[(X->SeqLen-1)%3] *3) / Width), 5);
      
      PlotLine(((i == 0) ? 0 : i-1),i, -3, -3,
               (((double)pGC[(X->SeqLen-2)%3]*3) / pWidth),
               (((double)GC[(X->SeqLen-2)%3] *3) / Width), 5);
    }

    if (plotAT) {
        PlotLine(((i == 0) ? 0 : i-1),i, +4, +4,
	       ((double)(pT_AT) / (pA_AT+pT_AT)),
	       ((double)(T_AT)  / (A_AT+T_AT)), 9);
	PlotLine(((i == 0) ? 0 : i-1),i, -4, -4,
               ((double)(pA_AT) / (pA_AT+pT_AT)),
	       ((double)(A_AT)  / (A_AT+T_AT)), 9);
    }
  }
}

// ------------------
//  Post analyse.
// ------------------
void SensorPlotter :: PostAnalyse(Prediction *pred)
{
}
