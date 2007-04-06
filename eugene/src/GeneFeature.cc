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
#include "GeneFeatureSet.h"
// -----------------------
//  Default constructeur
// -----------------------
GeneFeature::GeneFeature()
{
  id_=".";
  seqid_=".";
  source_=".";
  type_=".";
  locus_=NULL;
  score_=0.0;
  phase_='.';
  attributes_=NULL;
  valid_= true ;
}


//Constructor from a gff3 line

GeneFeature::GeneFeature ( char * line) 
{
  valid_= true ;
  line_=to_string(line);
  StringUtils::Chomp( line_ );
  strcpy (line , line_.c_str());
  if ( line_.at(0) != '#' )
  {  
    ParseLine(line); 
  }
}

// ---------------------------------------
//  ParseLine : Parse a gff3 line with a token, check SO/SOFA code and Parent Attribut.
// ---------------------------------------

void GeneFeature::ParseLine ( char * line ) 
{
  char * token = strtok (line,"\t");
  int i =0;
  int start, end;
  char strand;
  attributes_=NULL;
  locus_=NULL;
  while ( token != NULL )
  {
    //cout << "Token "<<i<< " valeur: "<<token<<endl;
    switch (i)
    {
      case 0 : 
      { //sequence id
	seqid_= to_string (token);
	break;
      }
      case 1 : 
      { //source
	source_= to_string (token);
	break;
      }
      case 2 : 
      { //type : SOFA/SO
	type_= to_string (token);
	if (  ( ! GeneFeatureSet::soTerms_->existsId(type_) ) && ( !GeneFeatureSet::soTerms_->existsName(type_) ) )  
	{
	  cout << "WARNING : " << type_ << " is not SO/SOFA referenced terms"<<endl;
	  valid_=false;
	}
	break;
      }
      case 3 : 
      {
	start= atoi (token);
	break;
      }
      case 4 : 
      {
	end= atoi (token);
	break;
      }
      case 5 : 
      {
	score_= atof (token);
	break;
      }
      case 6 : 
      {
	strand= token[0];
	break;
      }
      case 7 : 
      {
	phase_= token[0];
	break;
      }
      case 8 : 
      {
	//cout << "Attributes : " << token <<endl;
	attributes_ = new Attributes (to_string (token));
	id_=attributes_->getId();
	if ( id_ == "" )
	{
	  cout <<"WARNING : Feature has no ID >" << id_ << "<"<<endl;
	  valid_=false;
	}
	break;
      }
    }
    token = strtok (NULL,"\t");
    i++;
  }
  locus_  = new Locus (start, end, strand);
}

// -----------------------
//  Destructeur
// -----------------------
GeneFeature::~GeneFeature ( ) 
{ 
  delete locus_;
}

// -----------------------
// Accessor
// -----------------------

// Sequence Id
//
void GeneFeature::setSeqId ( string seqid ) 
{
  seqid_ = seqid;
}
string GeneFeature::getSeqId ( )
{
  return seqid_;
}

// Id of feature , copy from attributes
//

string GeneFeature::getId ( )
{
  return id_;
}

// Parent attributes get from attributes class if exists
//
string GeneFeature::getParent()
{
  string res="";
  if (attributes_ != NULL )
  {
    res= attributes_->getParent();
  }
  return res;
}

// Type
//
string GeneFeature::getType()
{
  return type_;
}

float GeneFeature::getScore()
{
  return score_;
}
// Valid attribut.
//
void GeneFeature::setValid (bool valid)
{
  valid_= valid;
}
bool GeneFeature::getValid ( )
{
  return valid_;
}
// ---------------------------------------
//  getString : Build a string gff3 format
// ---------------------------------------
string GeneFeature::getString ( )
{
  string geneFeature =seqid_+"\t"+source_+"\t"+type_+"\t";
  if ( locus_ != NULL )
  {
    geneFeature += to_string(locus_->getStart()) + "\t" + to_string(locus_->getEnd()) + "\t" + to_string(score_) + "\t" + to_string(locus_->getStrand()) + "\t";
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

//Locus
Locus * GeneFeature::getLocus ( )
{
  return locus_;
}