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
#include "Attributes.h"
#include "GeneFeatureSet.h";
//Attribut de classe;
    std::vector<std::string> Attributes::attributeNames_;
    
// -----------------------
//  Default constructeur
// -----------------------

Attributes::Attributes()
{
  InitializeDefault();
}



Attributes::Attributes ( std::string line ) 
{

  if (attributeNames_.empty()) {
    InitializeAttributeNames();
  }	
  //example : ID=match00001;Note=0907A18;Target=sp_P54263_SYN_THETH 6 34;
 
  InitializeDefault();
  char * oneAttributeString;
  char att[MAX_LINE];
  strcpy(att,line.c_str());

  oneAttributeString = strtok (att,"=;");
    char targetString[MAX_LINE];
    strcpy (targetString,"");
    int isFullLength=-1;
    int targetLength=-1;
    string targetSequence="";
    // Parse attributes 
    //  
    while ( oneAttributeString != NULL )
    {
      //Note=0907A18;
      char value[MAX_LINE]="";
      int i=0;
      while (i<NB_ATTRIBUTES_NAMES_ && strcmp (value,"") ==0 && oneAttributeString != NULL)  {
	
	if ( strcmp (attributeNames_[i].c_str(),oneAttributeString) ==0  ) 
	{
	  oneAttributeString = strtok (NULL,"=;");
	  strcpy(value,oneAttributeString);
	  // else at end of loop, add 1, wich no reason.
	  i--;
        }
	i++;
      } //end while
      if  ( i >= NB_ATTRIBUTES_NAMES_){
        break;
      }
      if ( attributeNames_[i] == "ID" && ( strcmp(value,"")!=0 ) ) {
	id_=to_string(value);
      }
      if ( attributeNames_[i] == "Name" && ( strcmp(value,"")!=0 )){
	name_=to_string(value);
      }
      if ( attributeNames_[i] == "Alias" && ( strcmp(value,"")!=0 )){
	alias_=to_string(value);
      }
      if ( attributeNames_[i] == "Parent" && ( strcmp(value,"")!=0 )){
	parent_=to_string(value);
	
      }
      if ( attributeNames_[i] == "Target" && ( strcmp(value,"")!=0 )) {
	strcpy (targetString,value);
      }
      if ( attributeNames_[i] == "Gap" && ( strcmp(value,"")!=0 )) {
	//gap_=value;
      }
      if ( attributeNames_[i] == "Derives_from" && ( strcmp(value,"")!=0 )) {
	derivesFrom_=to_string(value);
      }
      if ( attributeNames_[i] == "Note" && ( strcmp(value,"")!=0 )){
	note_=to_string(value);
      }
      if ( attributeNames_[i] == "Dbxref" && ( strcmp(value,"")!=0 )){
	dbxref_=to_string(value);
      }
      if ( attributeNames_[i] == "Ontology_term" && ( strcmp(value,"")!=0 )) {
	ontologyTerm_=to_string(value);
      }
      if ( attributeNames_[i] == "database" && ( strcmp(value,"")!=0 )) {
	database_=to_string(value);
      }
      if ( attributeNames_[i] == "is_full_length" && ( strcmp(value,"")!=0 )) {
	isFullLength=atoi(value);
      }
      if ( attributeNames_[i] == "target_length" && ( strcmp(value,"")!=0 )) {
	targetLength=atoi(value);
      }
      if ( attributeNames_[i] == "target_sequence" && ( strcmp(value,"")!=0 )) {
	targetSequence=to_string(value);
      }
      if ( attributeNames_[i] == "nb_exon" && ( strcmp(value,"")!=0 )) {
	nbExon_= atoi(value);
      }
      if ( attributeNames_[i] == "length" && ( strcmp(value,"")!=0 )) {
	length_= atoi(value);
      }
      oneAttributeString = strtok (NULL,"=;");
    } //end for all substring
    
    // check target attributes:
    if ( strcmp(targetString,"")!=0 ) { 
      char targetName[MAX_LINE];
      char targetStart[MAX_LINE],targetEnd[MAX_LINE],targetStrand=' ';
      //name
      char* targetAtt = strtok (targetString," ");
      strcpy(targetName,targetAtt);
      //start
      targetAtt = strtok (NULL," ");
      strcpy(targetStart,targetAtt);
      //end
      targetAtt = strtok (NULL," ");
      strcpy(targetEnd,targetAtt);
      //end
      targetAtt = strtok (NULL," ");
      if (targetAtt != NULL)
      {
	targetStrand=targetAtt[0];
      }
      target_ = new Target (targetName, targetSequence, targetLength, isFullLength, atoi(targetStart), atoi(targetEnd), targetStrand);
    }
}

// -----------------------
//  Destructeur
// -----------------------
Attributes::~Attributes ( ) 
{ 
// gaps_.clear();
}
//Initialize
void Attributes::InitializeDefault ( void)
{
  id_="";
  name_="";
  alias_="";
  note_="";
  parent_="";
  target_=NULL;
  database_="";
  derivesFrom_="";
  ontologyTerm_="";
  length_=-1;
  nbExon_=-1;
  //gaps_ = new std::vector();
}

// -----------------------
//Accessor
// -----------------------

// Id
//
void Attributes::setId ( std::string id ) 
{
	id_ = id;
}
std::string Attributes::getId ( )
{
	return id_;
}

std::string Attributes::getParent ( )
{
  return parent_;
}

std::string Attributes::getOntologyTerm ()
{
  return ontologyTerm_;
}

//getString
std::string Attributes::getString ( )
{	
	std::string my_attributes ="";
	if ( id_ != "" ) {
		my_attributes +=attributeNames_[0]+"="+id_+";";
	}
	if ( name_ != "" ) {
	  my_attributes +=attributeNames_[1]+"="+name_+";";
	}
	if ( alias_ != "" ) {
	  my_attributes +=attributeNames_[2]+"="+alias_+";";
	}
 	if ( parent_ != ""  ) {
	  my_attributes +=attributeNames_[3]+"="+parent_+";";
 	}
	if ( target_ != NULL ) {
		my_attributes +=target_->getString();
	}
/*	if ( !gaps_.empty() ) {
		//my_attributes +=gaps_.getString();
	}*/
	if ( derivesFrom_ != ""  ) {
	  my_attributes +=attributeNames_[6]+"="+derivesFrom_+";";
	}
	if ( note_ != "" ) {
	  my_attributes +=attributeNames_[7]+"="+note_+";";
	}
	if ( dbxref_ != "" ) {
	  my_attributes +=attributeNames_[8]+"="+dbxref_+";";
	}
	if ( ontologyTerm_ != "" ) {
	  my_attributes +=attributeNames_[9]+"="+ontologyTerm_+";";
	}
	if ( database_ != "" ) {
	  my_attributes +=attributeNames_[10]+"="+database_+";";
	}
	if ( nbExon_ != -1 ) {
	  my_attributes +=attributeNames_[14]+"="+to_string(nbExon_)+";";
	}
	if ( length_ != -1 ) {
	  my_attributes +=attributeNames_[15]+"="+to_string(length_)+";";
	}
	return my_attributes;
	
} 
