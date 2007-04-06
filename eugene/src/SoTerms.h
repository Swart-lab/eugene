
#ifndef SOTERMS_H
#define SOTERMS_H

#include <string>
#include <map>
#include "Const.h"
#include "System.h"
#include <iostream>
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
    SoTerms();
    virtual ~SoTerms ( );
    void loadFile( char * filename );
    bool existsId   (string name);
    bool existsName (string name);
    int size();
    void Print ();
    string getIdFromName (string name);

};
#endif // SOTERMS_H
