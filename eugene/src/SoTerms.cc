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
#include "SoTerms.h"
#include "Const.h"
// -----------------------
//  Default constructeur
// -----------------------
SoTerms::SoTerms( char * filename )
{
  FILE *fp;
  if (!(fp = fopen(filename, "r"))) {
    fprintf(stderr, "Class SoTerms : cannot open start file %s\n", filename);
    exit(2);
  }
  char tag[60] , value[100];
  
  int i=1;
  int j=0;
  char  line[MAX_LINE];
  while(fp  &&  fgets (line, MAX_LINE, fp) != NULL) {
    j++;
    if (line[0] == 'i' && line[1] == 'd') {
      i = sscanf(line, "id: %s", &value);
      if (i > 0) {
	char soId[60];
	strcpy (soId, value );
	fgets (line, MAX_LINE, fp);
	i = sscanf(line, "name: %s", &value);
	idToName_[soId]=value;
	nameToId_[value]=soId;
      }
    }
  }
  fclose(fp);
}

SoTerms::~SoTerms ( )
{}

void SoTerms::Print (void)
{
  map<string,string>::iterator it;
  for ( it=idToName_.begin() ; it != idToName_.end(); it++ )
  {
    cout << "idToName_ Key : " << (*it).first << " Value : " << (*it).second << endl;
  }
}

bool SoTerms::existsId(string id)
{
  map<string,string>::iterator it=idToName_.find(id);
  if (it == idToName_.end()) {
      return false;
  }
  else {
    return true;
  }
}

bool SoTerms::existsName( string name)
{
  map<string,string>::iterator it=nameToId_.find(name);
  if (it == idToName_.end()) {
    return false;
  }
  else {
    return true;
  }
}