// ------------------------------------------------------------------
// Copyright (C) 2009 INRA <eugene@ossau.toulouse.inra.fr>
//
// This program is open source; you can redistribute it and/or modify
// it under the terms of the Artistic License (see LICENSE file).
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
//
// You should have received a copy of Artistic License along with
// this program; if not, please see http://www.opensource.org
//
// $Id
// ------------------------------------------------------------------
// File:     State.h
// Contents: Class State
// ------------------------------------------------------------------

#ifndef  STATE_H_INCLUDED
#define  STATE_H_INCLUDED

#include "Prediction_cte.h"

class State
{
private:
    char state;
public:
    State();
    State(char state);
    char GetState();
    short int GetFrame(); /* Return the frame - See the definition on the EuGene Trac */
    bool IsIntergenic(void);
    bool IsIntron(void);
    bool IsForwardIntron(void);
    bool IsReverseIntron(void);
    bool IsUTR(void);
    bool IsUTR5(void);
    bool IsUTR3(void);
    bool IsCodingExon(void);
    bool IsForwardCodingExon(); /* coding exon on the strand forward */
    bool IsReverseCodingExon(); /* coding exon on the strand reverse */
    bool IsTranscribedAndUnspliced(void);
    bool InStartStopRegion(); /* true if its an element including between a start and a stop codon */
    bool IsDefined();

};

#endif
