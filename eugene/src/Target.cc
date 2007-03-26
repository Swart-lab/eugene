/***************************************************************************
 *   Copyright (C) 2007 by Noirot   *
 *   cnoirot@milieu   *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/
#include "Target.h"

// -----------------------
//  Default constructeur
// -----------------------
Target::Target ( ) 
{
  name_ = "";
  range_ = NULL;
  sequenceData_ = "";
  targetLength_ = 0;
  isFullLength_ = UNKNOWN; 

}

Target::Target (std::string name, Range * range , std::string sequenceData,  int targetLength, int isFullLength) 
{
  name_ = name;
  range_ = new Range (range);
  sequenceData_ = sequenceData;
  targetLength_ = targetLength;
  isFullLength_ = isFullLength; 
}

// -----------------------
//  Destructeur
// -----------------------

Target::~Target ( ) 
{ 
	delete range_;
}


// -----------------------
//Accessor
// -----------------------

// Name
//
void Target::setName ( std::string name ) 
{
 	name_ = name;
}
std::string Target::getName ( )
{
	return name_;
}

// Range
//
void Target::setRange ( Range * range ) 
{
 	range_ = range;
}
Range * Target::getRange ( )
{
	return range_;
}

// SequenceData_
//
void Target::setSequenceData ( std::string sequenceData ) 
{
 	sequenceData_ = sequenceData;
}
std::string Target::getSequenceData ( )
{
	return sequenceData_;
}


//targetLength_
//
void Target::setGlobalLength ( int targetLength )
{
	targetLength_=targetLength;
}
int Target::getGlobalLength ( )
{
	return targetLength_;
}


//isFullLength_
//
void Target::setIsFullLength ( int isFullLength )
{
	isFullLength_=isFullLength;
}
int Target::getIsFullLength ( )
{
	return isFullLength_;
}


//getString
std::string Target::getString ( )
{
	std::string my_target = "";
	if ( name_ != ""  && range_ != NULL)
	{	
		 my_target = "target="+name_+" "+range_->getString()+";";
	}
	if ( isFullLength_ != UNKNOWN )
	{
		my_target += "is_full_length=" + to_string(isFullLength_) + ";";
	}
	if ( targetLength_ != 0 )
	{
		my_target += "target_length=" + to_string(targetLength_)+ ";";
	}
	if ( sequenceData_ != ""  )
	{
		my_target += "target_sequence=" + sequenceData_+ ";";
	}
	return my_target;
}
