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

// -----------------------
//  Default constructeur
// -----------------------
Attributes::Attributes()
{
  id_="";
  name_="";
  alias_="";
  note_="";
  //parent_=NULL;
  target_=NULL;
  database_="";
  //deriveFrom_=NULL;
  ontologyTerm_="";
  length_=-1;
  nbExon_=-1;
  //gaps_ = new std::vector();
}

Attributes::Attributes ( std::string line ) 
{
	if (attributeNames_.empty()) {
		InitAttributeNames();
	}	
	std::vector<std::string> subStrings;
	//example : ID=match00001;Note=0907A18;Target=sp_P54263_SYN_THETH 6 34;
	int nbAttribute = StringUtils::SplitString(line ,";",subStrings);
	if(subStrings.empty()) { 
		Attributes();
	}
	else {
		std::vector<std::string>::const_iterator iSubStrings; // Construction d'un itérateur pour parcourir les éléments du vecteur
		string targetString="";
		int isFullLength=-1;
		int targetLength=0;
		string targetSequence="";
		
		// Parse attributes 
		for(iSubStrings=subStrings.begin();iSubStrings!=subStrings.end();iSubStrings++)
		{
			int pos=-1;
			string value="";
			int i=0;
			while (i<NB_ATTRIBUTES_NAMES_ && value=="")
			{
				if ( (pos=(*iSubStrings).find (attributeNames_[i],0)) != string::npos  ) {
					value=(*iSubStrings).substr(pos+attributeNames_[i].length(), (*iSubStrings).length());
					cout << (*iSubStrings)<< "att name " << attributeNames_[i] << "i:" <<i <<endl;
					// else at end of loop, add 1, wich no reason.
					i--;
				}
 				i++;
			}

			if ( attributeNames_[i] == "ID=" && value !="") {
				id_=value;
			}
			if ( attributeNames_[i] == "Name=" && value !=""){
				name_=value;
			}
			if ( attributeNames_[i] == "Alias=" && value !=""){
				alias_=value;
			}
			if ( attributeNames_[i] == "Parent=" && value !=""){
				//alias_=value;
			}
			if ( attributeNames_[i] == "Target=" && value !="") {
				targetString=value;
			}
			if ( attributeNames_[i] == "Gap=" && value !="") {
				//gap_=value;
			}
			if ( attributeNames_[i] == "Derives_from=" && value !="") {
				//derives_from_=value;
			}
			if ( attributeNames_[i] == "Note=" && value !=""){
				note_=value;
 			}
			if ( attributeNames_[i] == "Dbxref=" && value !=""){
				dbxref_=value;
			}
			if ( attributeNames_[i] == "Ontology_term=" && value !="") {
				ontologyTerm_=value;
			}
			if ( attributeNames_[i] == "database=" && value !="") {
				database_=value;
			}
			if ( attributeNames_[i] == "is_full_length=" && value !="") {
				isFullLength=atoi(value.c_str());
			}
			if ( attributeNames_[i] == "target_length=" && value !="") {
				targetLength=atoi(value.c_str());
			}
			if ( attributeNames_[i] == "target_sequence=" && value !="") {
				targetSequence=value;
			}
			if ( attributeNames_[i] == "nb_exon=" && value !="") {
				nbExon_= atoi(value.c_str());
			}
			if ( attributeNames_[i] == "length=" && value !="") {
				length_= atoi(value.c_str());
			}
		} 
		
		// check target attributes:
		if ( targetString!="") { 
			//isFullLength targetLength targetSequence
			std::vector<std::string> targetAtt;
			int targetNbAtt = StringUtils::SplitString(targetString ," ",targetAtt);
			cout << "Target "<< targetString << " nb el: " << targetNbAtt <<endl;
			char strand=' ';
			if (targetNbAtt == 3 ) {
				strand=(targetAtt[3].c_str())[0];
			}
			//Target ( std::string name, std::string sequenceData,  int targetLength, int isFullLength, int start, int end, char strand);
			target_ = new Target (targetAtt[0],targetSequence,targetLength,isFullLength,atoi(targetAtt[1].c_str()), atoi(targetAtt[2].c_str()), strand);
		}
		
	}

}

// -----------------------
//  Destructeur
// -----------------------
Attributes::~Attributes ( ) 
{ 
 gaps_.clear();
}


//getString
std::string Attributes::getString ( )
{
	std::string my_attributes ="";
	if ( id_ != "" ) {
		my_attributes +=attributeNames_[0]+id_+";";
	}
	if ( name_ != "" ) {
		my_attributes +=attributeNames_[1]+name_+";";
	}
	if ( alias_ != "" ) {
		my_attributes +=attributeNames_[2]+alias_+";";
	}
// 	if ( parent_ != NULL ) {
// 		my_attributes +=attributeNames_[3]+parent_.getString()+";";
// 	}
	if ( target_ != NULL ) {
		my_attributes +=target_->getString();
	}
	if ( !gaps_.empty() ) {
		//my_attributes +=gaps_.getString();
	}
// 	if ( derivesFrom_ != NULL ) {
// 		my_attributes +=attributeNames_[6]+derivesFrom_.getId()+";";
// 	}
	if ( note_ != "" ) {
		my_attributes +=attributeNames_[7]+note_+";";
	}
	if ( dbxref_ != "" ) {
		my_attributes +=attributeNames_[8]+dbxref_+";";
	}
	if ( ontologyTerm_ != "" ) {
		my_attributes +=attributeNames_[9]+ontologyTerm_+";";
	}
	if ( database_ != "" ) {
		my_attributes +=attributeNames_[10]+database_+";";
	}
	if ( length_ != -1 ) {
		my_attributes +=attributeNames_[11]+to_string(length_)+";";
	}
	if ( nbExon_ != -1 ) {
		my_attributes +=attributeNames_[12]+to_string(nbExon_)+";";
	}
	return my_attributes;
} 