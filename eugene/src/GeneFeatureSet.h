
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
    vector <GeneFeature *> vRefFeatures_;
    map <string, int> mPosFeatures_;
    int lastIndex_;
    map <string, vector<string> > parentToChildren_;
  public:
    static SoTerms * soTerms_ ;
    // Constructors/Destructors
    GeneFeatureSet ( );
    GeneFeatureSet ( char* featuresFilename, char* soTermsFilename );
    GeneFeatureSet ( char* soTermsFilename );  
    virtual ~GeneFeatureSet ( );
    bool existsGeneFeature ( string geneFeatureId ) ;
    void printFeature();
    vector< GeneFeature *>::iterator getIterator ();
    map<string, vector<string> >::iterator getIteratorParentToChildren ();
    int getNbFeature ();
    int getNbParentFeature ();
    GeneFeature * getGeneFeature (string id );
};

#endif // GENEFEATURESET_H
