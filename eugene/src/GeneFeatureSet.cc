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
  lastIndex_=0;
  // si premier object geneFeatureSet il faut initialiser la variable de classe
  if ( GeneFeatureSet::soTerms_->size() == 0 ) {
    //cout<<"Load SoTerms"<<endl;
    GeneFeatureSet::soTerms_->loadFile(soTermsFilename);
  }
  FILE *fp;
  if (!(fp = fopen(featuresFilename, "r"))) {
    fprintf(stderr, "Class GeneFeatureSet : cannot open file %s\n", featuresFilename);
  }
  else
  {
  char value[100]="";
  int i=1;
  char line[MAX_LINE];
  while(fp  &&  fgets (line, 1500, fp) != NULL) 
  {
    if ( line[0] != '#' )
    {  
      GeneFeature * tempGeneFeature = new GeneFeature (line);
      // cout << "Test >"<< tempGeneFeature->getId() << "< with parent >" << tempGeneFeature->getParent()<<"<" <<endl;

      if ( tempGeneFeature->getParent()!="" && ! existsGeneFeature( tempGeneFeature->getParent() ) )
      {
	tempGeneFeature->setValid(false);
	cout << "WARNING : parent >" << tempGeneFeature->getParent()  << "< does not exist "<<endl;
      }
      else
      {
	if (tempGeneFeature->getParent()!="")
	{
	   parentToChildren_[tempGeneFeature->getParent()].push_back(tempGeneFeature->getId());
	}
      }
     
      if ( tempGeneFeature->getValid() )
      {
	vRefFeatures_.push_back( tempGeneFeature);
	mPosFeatures_[tempGeneFeature->getId()]=lastIndex_;
	lastIndex_++;
      }
     
    }
  }
  fclose(fp);
  }

}

//Constructor
GeneFeatureSet::GeneFeatureSet (char* soTermsFilename ) 
{
 
  // si premier object geneFeatureSet il faut initialiser la variable de classe
  if ( GeneFeatureSet::soTerms_->size() == 0 ) {
    //cout<<"Load SoTerms"<<endl;
    GeneFeatureSet::soTerms_->loadFile(soTermsFilename);
  }

}


/*void GeneFeatureSet::getIteratorParentToChildren()
{
  map <string, GeneFeature *>::iterator it;
  map <string, vector > analyse;
  for ( it=features_.begin() ; it != features_.end(); it++ )
  {
    cout << "features_ Key : >" << (*it).first << "< Value : >" << *((*it).second)<<"<";
    if ( ((*it).second)->getParent() != "" )
    {
      analyse[((*it).second)->getParent()]
    }
  }
 // cout << "features_ Key : " << (*it).first << " Value : " << (*it).second << endl;
}
*/
// -----------------------
//  Destructeur
// -----------------------
GeneFeatureSet::~GeneFeatureSet ( ) 
{ 
  cout << "detroy GeneFeatureSet"<<endl;
  vector < GeneFeature *>::iterator it;
  for ( it=vRefFeatures_.begin() ; it != vRefFeatures_.end(); it++ )
  {
    delete *it ;
  }
  vRefFeatures_.clear();
  mPosFeatures_.clear();
}

bool GeneFeatureSet::existsGeneFeature ( string geneFeatureId ) 
{
  map<string, int>::iterator it=mPosFeatures_.find(geneFeatureId);
  bool res = true;
  if ( mPosFeatures_.end() == it )
    res=false;
  return res;
}

map<string, vector<string> >::iterator GeneFeatureSet::getIteratorParentToChildren ()
{
  return parentToChildren_.begin() ;
}

int GeneFeatureSet::getNbParentFeature ()
{
  return parentToChildren_.size();
}

vector< GeneFeature *>::iterator GeneFeatureSet::getIterator ()
{
  return vRefFeatures_.begin() ;
}

int GeneFeatureSet::getNbFeature ()
{
  return vRefFeatures_.size();
}

GeneFeature * GeneFeatureSet::getGeneFeature (string id )
{
  return vRefFeatures_[(*mPosFeatures_.find(id)).second];
}

void GeneFeatureSet::printFeature()
{
  int i = 0 ;
  for ( i=0 ; i < lastIndex_ ; i++ )
  {
    cout << *(vRefFeatures_[i]) << endl;
  }
}
