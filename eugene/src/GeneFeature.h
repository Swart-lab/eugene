
#ifndef GENEFEATURE_H
#define GENEFEATURE_H

#include <string>
#include <iostream>

#include "StringUtils.h"

//Feature class:
#include <Range.h>
#include <Attributes.h>
#include <GeneFeatureSet.h>
using namespace std; 

/**
 * class GeneFeature
 */


class GeneFeature
{
  private:
    string id_;
    string seqid_;
    string source_;
    string type_;
    Range * range_;
    float score_;
    char phase_;
    Attributes * attributes_;
    
  public:
  
    string line_;
    // Constructors/Destructors
    GeneFeature ( );
    GeneFeature ( string line ); 
    virtual ~GeneFeature ( );
    void setSeqId ( string id );
    string getSeqId ( );
    string getString ( );
    
    // Public attribute accessor methods
    //  
    friend ostream& operator<<(ostream& o, GeneFeature geneFeature )
    {
      return o << geneFeature.getString()<< endl;
    }

};

#endif // GENEFEATURE_H
