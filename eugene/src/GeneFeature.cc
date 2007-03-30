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
#include "GeneFeature.h"

// -----------------------
//  Default constructeur
// -----------------------
GeneFeature::GeneFeature()
{
  id_=".";
  seqid_=".";
  source_=".";
  type_=".";
  range_=NULL;
  score_=0.0;
  phase_='.';
  attributes_=NULL;
}

GeneFeature::GeneFeature ( string line ) 
{
  line_=line;
  vector<string> subStrings;
  
  //example : ID=match00001;Note=0907A18;Target=sp_P54263_SYN_THETH 6 34;
  int nbColumn = StringUtils::SplitString(line ,"\t",subStrings);
  cout << "nb column " << nbColumn << endl;
  if(subStrings.empty() || nbColumn < 7 || nbColumn > 9) { 
    cout << "Error : nb column " << nbColumn << endl;
    GeneFeature();
  }
  else {
    //check sofa / so terms
    if (  ! GeneFeatureSet::soTerms_->existsId(subStrings[2]) && !GeneFeatureSet::soTerms->existsName(subStrings[2]) )  {
      cout << "WARNING : " << subStrings[2] << " is not SO/SOFA referenced terms"<<endl;
    }
    // Parse attributes 
    seqid_  = subStrings[0];
    source_ = subStrings[1];
    type_   = subStrings[2];
    score_  = atof(subStrings[5].c_str());
    range_  = new Range (atoi(subStrings[3].c_str()), atoi(subStrings[4].c_str()),(subStrings[6].c_str())[0]);
    phase_  = (subStrings[7].c_str())[0];
    if ( nbColumn == 8 && subStrings[8] != "") {
      attributes_ = new Attributes (subStrings[8]);
      id_=attributes_->getId();
    }
    else {
      attributes_ = NULL;
    }
  } 
}

// -----------------------
//  Destructeur
// -----------------------
GeneFeature::~GeneFeature ( ) 
{ 

}

// -----------------------
//Accessor
// -----------------------

// Id
//
void GeneFeature::setSeqId ( string seqid ) 
{
  seqid_ = seqid;
}
string GeneFeature::getSeqId ( )
{
  return seqid_;
}



//getString
string GeneFeature::getString ( )
{
  string geneFeature =seqid_+"\t"+source_+"\t"+type_+"\t";
  if ( range_ != NULL )
  {
    geneFeature += to_string(range_->getStart()) + "\t" + to_string(range_->getEnd()) + "\t" + to_string(score_) + "\t" + to_string(range_->getStrand()) + "\t";
  }
  else
  {
    geneFeature += ".\t.\t" + to_string(score_) + "\t.\t";
  }
  
  geneFeature += to_string(phase_)+"\t";
  if ( attributes_ != NULL )
  {
    geneFeature +=attributes_->getString();
  }
  return geneFeature;
} 
