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
#include "GeneFeatureSet.h"

// -----------------------
//  Default constructeur
// -----------------------
GeneFeatureSet::GeneFeatureSet()
{}

//Static attribut
SoTerms *  GeneFeatureSet::soTerms_= new SoTerms ();

//Constructor
GeneFeatureSet::GeneFeatureSet ( char* featuresFilename, char* soTermsFilename ) 
{
  // si premier object geneFeatureSet il faut initialiser la variable de classe
  if ( GeneFeatureSet::soTerms_->size() == 0 ) {
    //cout<<"Load SoTerms"<<endl;
    GeneFeatureSet::soTerms_->loadFile(soTermsFilename);
  }

  FILE *fp;
  if (!(fp = fopen(featuresFilename, "r"))) {
    fprintf(stderr, "Class GeneFeatureSet : cannot open file %s\n", featuresFilename);
    exit(2);
  }
  char value[100]="";
  int i=1;
  char line[MAX_LINE];
  while(fp  &&  fgets (line, 1500, fp) != NULL) {
    //cout  << "LINE: "<< line;
    if ( line[0] != '#' )
    {  
      
      GeneFeature * tempGeneFeature = new GeneFeature (line);
      if ( tempGeneFeature->getParent()!="" && ! existsGeneFeature( tempGeneFeature->getParent() ) )
      {
	tempGeneFeature->setValid(false);
	cout << "WARNING : parent " << tempGeneFeature->getParent()  << " does not exist "<<endl;
      }
      if ( tempGeneFeature->getValid() )
      {
	features_[tempGeneFeature->getId()]=tempGeneFeature;
	//cout << "KEEP : " << tempGeneFeature->getId()<<  endl;
      }
     
    }
    
  }
  fclose(fp);
  //cout<<"End constructor : "<<featuresFilename <<" SO term : "<< soTermsFilename <<endl;
}


void GeneFeatureSet::printFeature()
{
  map <string, GeneFeature *>::iterator it;
  for ( it=features_.begin() ; it != features_.end(); it++ )
  {
    cout << "features_ Key : " << (*it).first << " Value : " << *((*it).second);
    
  }
 // cout << "features_ Key : " << (*it).first << " Value : " << (*it).second << endl;
}


// -----------------------
//  Destructeur
// -----------------------
GeneFeatureSet::~GeneFeatureSet ( ) 
{ 
    features_.clear();
}

bool GeneFeatureSet::existsGeneFeature ( string geneFeatureId ) 
{
  map<string, GeneFeature *>::iterator it=features_.find(geneFeatureId);
  bool res = true;
  if (it == features_.end()) {
    res=false;
  }
  return res;
}

map<string, GeneFeature *>::iterator GeneFeatureSet::getIterator ()
{
  return features_.begin() ;
}

int GeneFeatureSet::getNbFeature ()
{
  return features_.size();
}