
#ifndef SOTERMS_H
#define SOTERMS_H

#include <string>
#include <map>
#include "StringUtils.h"
using namespace std; 

/**
 * class SoTerms
    
    example of enregistrement
    [Term]
    id: SO:0000000
    name: Sequence_Ontology
    subset: SOFA

 */


class SoTerms
{
  private:
    map <string,string> idToName_;
    map <string,string> nameToId_;
    
  public:
  
    // Constructors/Destructors
    //SoTerms ( );
    SoTerms( char * filename );
    virtual ~SoTerms ( );
    bool existsId   (string name);
    bool existsName (string name);
    void Print ();

};
#endif // SOTERMS_H
