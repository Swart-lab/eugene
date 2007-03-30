
#ifndef ATTRIBUTES_H
#define ATTRIBUTES_H

#include <string>
#include <vector>
#include <iostream>

#include "StringUtils.h"

//Feature class:
#include <Target.h>
#include <Gap.h>
/**
  * class Attributes
  */


class Attributes
{


private:
 
  //From gff3 definition
  std::string id_;
  std::string name_;
  std::string alias_;
  std::string parent_;
  Target * target_;
  // on an alignment list of the Gap.
  std::vector<Gap> gaps_;
  std::string derivesFrom_;
  std::string note_;
  std::string dbxref_;
  std::string ontologyTerm_;
  std::string database_;  
  
  //For eugene definition
  int length_;
  int nbExon_;

public:
  

  // Constructors/Destructors
  Attributes ( );
  Attributes ( std::string line ); 
  void Attributes::InitializeDefault ( void);
  virtual ~Attributes ( );
  void setId ( std::string id );
  std::string getId ( );
  std::string getString ( );

  // Puplic attribute accessor methods
  //  
  friend std::ostream& operator<<(std::ostream& o, Attributes attributes )
  {
	return o << attributes.getString();
  }
	

  //Should be static members , 
  static const int NB_ATTRIBUTES_NAMES_=16;
  std::vector<std::string> attributeNames_; // "ID=","Note=","Target="
  void InitializeAttributeNames()
  {
	attributeNames_.push_back("ID=");
	attributeNames_.push_back("Name=");
	attributeNames_.push_back("Alias=");
	attributeNames_.push_back("Parent=");
  	attributeNames_.push_back("Target=");
	attributeNames_.push_back("Gap=");
	attributeNames_.push_back("Derives_from=");
	attributeNames_.push_back("Note=");
	attributeNames_.push_back("Dbxref=");
	attributeNames_.push_back("Ontology_term=");
	attributeNames_.push_back("database=");
	attributeNames_.push_back("is_full_length=");
	attributeNames_.push_back("target_length=");
	attributeNames_.push_back("target_sequence=");
	attributeNames_.push_back("nb_exon=");
	attributeNames_.push_back("length=");
  }

};


  // Private attributes
	//From GFF3 specification 
		//ID=match00001;Note=0907A18;Target=sp_P54263_SYN_THETH 6 34;Ontology_term=SO:000012;
	//From sensor : 
		//database=cDNAMt;
		//In class Target : is_full_length=1;target_length=231;
	//Eugene output: 
		//nb_exon=7;length=3598
#endif // ATTRIBUTES_H
