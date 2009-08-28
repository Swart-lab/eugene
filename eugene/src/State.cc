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

#include "State.h"

/*************************************************************
 **                      State                              **
 *************************************************************/

// ----------------------
//  Default constructor.
// ----------------------
State::State()
{
    this->state = -1;
}

// ----------------------
//  Constructor.with a parameter: the state in char format
// ----------------------
State::State(char state)
{
    this->state = state;
}

// ----------------------
//  Return the char value of the state
// ----------------------
char State::GetState()
{
    return this->state;
}

// ----------------------
//  Return true if it is intergenic
// ----------------------
bool State::IsIntergenic(void)
{
    if (State2Status[this->state] == UNTRANSCRIBED)
    {
        return true;
    }
    return false;
}


// ----------------------
//  Return true if it an intron
// ----------------------
bool State::IsIntron(void)
{
    if (State2Status[this->state] == SPLICED_TRANSCRIBED)
    {
        return true;
    }
    return false;
}

// ----------------------
//  Return true if it is an UTR
// ----------------------
bool State::IsUTR(void)
{
    if (State2Status[this->state] == UNTRANSLATED)
    {
        return true;
    }
    return false;
}

// ----------------------
//  Return true if it is an UTR in 5'
// ----------------------
bool State::IsUTR5(void)
{
    if (this->state == UTR5F  ||  this->state == UTR5R)
    {
        return true;
    }
    return false;
}

// ----------------------
//  Return true if it is an UTR in 3'
// ----------------------
bool State::IsUTR3(void)
{
    if (this->state == UTR3F  ||  this->state == UTR3R)
    {
        return true;
    }
    return false;
}

// ----------------------
//  Return true if it is a coding exon
// ----------------------
bool State::IsCodingExon(void)
{
    if (State2Status[this->state] == TRANSLATED)
    {
        return true;
    }
    return false;
}

// ----------------------
//  Return true if coding exon on the strand forward
// NOTE  A modifier! Return (IsCodingExon() && IsForward())
// ----------------------
bool State::IsForwardCodingExon()
{
    if (1 <= State2Frame[this->state] && State2Frame[this->state] <= 3)
    {
        return true;
    }
    return false;
}

// ----------------------
//  Return true if coding exon on the reverse forward
// NOTE  A modifier! Return (IsCodingExon() && IsReverse())
// ----------------------
bool State::IsReverseCodingExon()
{
    if (-3 <= State2Frame[this->state] && State2Frame[this->state] <= -1)
    {
        return true;
    }
    return false;
}

// ----------------------
//  Return true if it is transcribed and unspliced: UTR or coding exon
// ----------------------
bool State::IsTranscribedAndUnspliced(void)
{
    if (this->IsUTR() || this->IsCodingExon())
    {
        return true;
    }
    return false;
}

// -----------------------
// Return true if its an intron in the forward strand
// NOTE: remplacer ce test par IsIntron && IsForward
// -----------------------
bool State::IsForwardIntron(void)
{
    if (State2Frame[this->state] == FrameIntronF)
    {
        return true;
    }
    return false;
}

// -----------------------
// Return true if its an intron in the reverse strand
// NOTE: remplacer ce test par IsIntron && IsReverse
// -----------------------
bool State::IsReverseIntron(void)
{
    if (State2Frame[this->state] == FrameIntronR)
    {
        return true;
    }
    return false;
}


// -----------------------
// Return true if its an element including between a start and a stop codon
// Not IG, UTR and UTR Intron
// -----------------------
bool State::InStartStopRegion(void)
{
    if (this->state < InterGen)
        return true;
    return false;
}

// ------------------------
//  Return the frame - See the definition on the EuGene trac
// ------------------------
short int State::GetFrame()
{
    return State2Frame[this->state];
}

// ------------------------
//  Return true if the state is defined
// NOTE: for the moment, test if the value of state is positive; this test can be improved
// ------------------------
bool State::IsDefined()
{
    if (this->state >= 0)
        return true;
    return false;
}
