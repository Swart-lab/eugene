
#ifndef GENEFEATURESET_H
#define GENEFEATURESET_H
//#define _VACPP_TR1

#include <string>
#include <iostream>

#include "GeneFeature.h"
#include "SoTerms.h"
#include "Const.h"

using namespace std;

/**
 * class GeneFeatureSet
 */

    
class GeneFeatureSet
{
  
  private:
    map <string, GeneFeature *> features_;
  public:
    static SoTerms * soTerms_ ;
    // Constructors/Destructors
    GeneFeatureSet ( );
    GeneFeatureSet ( char* featuresFilename, char* soTermsFilename ); 
    virtual ~GeneFeatureSet ( );
    bool existsGeneFeature ( string geneFeatureId ) ;
    void printFeature();
    map<string, GeneFeature *>::iterator getIterator ();
    int getNbFeature ();
 
};

#endif // GENEFEATURESET_H
