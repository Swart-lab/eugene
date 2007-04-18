
#ifndef GENEFEATURE_H
#define GENEFEATURE_H

#include <string>
#include <iostream>
#include "StringUtils.h"
//Feature class:
#include "Locus.h"
#include "Attributes.h"
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
    Locus * locus_;
    double score_;
    char phase_;
    Attributes * attributes_;
    bool valid_ ;
    void ParseLine( char * line );
    
  public:
  
    string line_;
    // Constructors/Destructors
    GeneFeature ( );
    GeneFeature ( char * ); 
    //GeneFeature (string seqId, string source, string type, int start, int end, double score, char strand, char phase, char * attributes);
    virtual ~GeneFeature ( );
    
    void setValid ( bool valid);
    bool getValid ( );
    void setSeqId ( string id );
    
    string getSeqId ( );
    string getString ( );
    string getId ();
    string getParent();
    string getType();
    Locus * getLocus();
    double getScore();
    // Public attribute accessor methods
    //  
    friend ostream& operator<<(ostream& o, GeneFeature geneFeature )
    {
      return o << geneFeature.getString()<< endl;
    }

};

#endif // GENEFEATURE_H
