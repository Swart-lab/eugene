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
#include "Locus.h"

// -----------------------
//  Default constructeur
// -----------------------
Locus::Locus ( ) 
{
  start_ =0;
  end_   =0;
  strand_=' ';
}

Locus::Locus (int start, int end, char strand) 
{
  start_ =start;
  end_   =end;
  strand_=strand;
  cout << "Constructeur Locus :  " << getString() <<endl;
}

Locus::Locus (Locus *locus )
{
  start_ =locus->getStart();
  end_   =locus->getEnd();
  strand_=locus->getStrand();
}

// -----------------------
//  Destructeur
// -----------------------
Locus::~Locus ( ) { }


//getString
std::string Locus::getString ( )
{
	std::string my_locus = to_string(start_)+" "+to_string(end_);
	if ( strand_ != ' ' )
	{
		my_locus +=" "+to_string(strand_);
	}
	return my_locus;
}

// surcharge <<
// ostream& operator<<(ostream &o, const Locus *locus )
// {
// 	return o << locus->getString();
// }

// -----------------------
//  Start
// -----------------------
void Locus::setStart ( int start )
{
	start_=start;
}
int Locus::getStart ( )
{
	return start_;
}

// -----------------------
//  End
// -----------------------
void Locus::setEnd ( int end )
{
	end_=end;
}
int Locus::getEnd ( )
{
	return end_;
}


// -----------------------
// Strand
// -----------------------
void Locus::setStrand ( char strand )
{
	strand_=strand;
}
char Locus::getStrand ( )
{
	return strand_;
}


