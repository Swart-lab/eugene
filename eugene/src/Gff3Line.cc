/**  \file    Gff3Line.cpp
 *  \author  S�bastien Letort
 *  \date    10 novembre 2006
 */

#include "Gff3Line.h"

///The default source name
static const std::string EUGENE = "EuGene";
///an integer considered as undefined
static const int INDEFINI = -9999;

//------------------------
///a terme gerer les code avec un tableau/liste ?
std::string
  Sofa::getName(unsigned int code_so)
{
  switch(code_so)
  {//apparait dans le gene
    case SOFA_GENE :  return "gene";  break;
    case SOFA_EXON :  return "exon";  break;
    case SOFA_INTRON :  return "intron";  break;
    case SO_UTR_INTRON :  return "intron_utr";  break;
    case SOFA_INTERGEN :  return "intergenic_region";  break;
   //aparait dans l'ARNm (episse)
    case SOFA_MRNA :  return "mRNA";  break;
//    case SOFA_POLY_A_SEQ :  return "polyA_sequence";  break;
//    case SOFA_POLY_A_SITE :  return "polyA_site";  break;
    case SOFA_CDS :  return "CDS";  break;
    case SOFA_5_UTR :  return "five_prime_UTR";  break;
    case SOFA_3_UTR :  return "three_prime_UTR";  break;
    //par default on renvoie le code au format so::000xx
    default :  return Sofa::codeToString(code_so);
  }
}

std::string
  Sofa::codeToString(const int code_so)
{
  // utiliser un flux de sortie pour creer la chaine
  std::ostringstream oss;
  oss << "SO:" << std::setfill('0') << std::setw(7) << code_so;
  // renvoyer une string
  return oss.str();
}
//------------------------
//les constructeurs sont places la pour avoir acces a la constante EUGENE
///constructeur pour une sequence donnee
Gff3Line::Gff3Line(std::string seq_id)
      : seq_id_(seq_id), source_(EUGENE), type_sofa_(INDEFINI),
        start_(INDEFINI), end_(INDEFINI), score_(INDEFINI),
        strand_('.'), phase_(INDEFINI), attr_(0)
{
//fprintf(stderr, "ctr2 gff3Line\n");// << std::endl;
}

//renvoie le type sofa le + approprie au type eugene
int getTypeSofa(int state_egn, bool coding, bool sofa)
{
  int ret_val;
    if(state_egn <= InitR3)
//        str = "five_prime_coding_exon";//"E.Init";
        return (sofa) ? SOFA_EXON : (coding ? SO_5_CODING_EXON : SO_5_EXON);
    if(state_egn <= SnglR3)
//        str = "single_exon";//"E.Sngl";
        return (sofa) ? SOFA_EXON : SO_SINGLE_EXON;
    if(state_egn <= IntrR3)
        //interior_exon et interior_coding_exon ne sont pas dans SOFA
      return (sofa) ? SOFA_EXON : SO_INTERIOR_EXON;  //"E.Intr";
//        str = "SO:0000201";  //"E.Intr";
    if(state_egn <= TermR3)
//        str = "three_prime_coding_exon";//"E.Term";
      return (sofa) ? SOFA_EXON : (coding ? SO_3_CODING_EXON : SO_3_EXON);
    if(state_egn <= IntronR2AG)
      return SOFA_INTRON;    //"Intron";
    if(state_egn == InterGen)
      return SOFA_INTERGEN;  //"InterG";
    if(state_egn == UTR5F || state_egn == UTR5R)
      return SOFA_5_UTR;  //"UTR5";
    if(state_egn == UTR3F || state_egn == UTR3R)
      return SOFA_3_UTR;  //"UTR3";
//    if(state_egn >= IntronU5F)
    if(state_egn <= IntronU3R)
      return SO_UTR_INTRON;    //"Intron";
//     si ce n'est rien de tout ca, je renvoies SOFA_EXON
//     a cause des attributs,
//     je dois pouvoir recuperer un SOFA_EXON avec un sofa a false
//     a changer quand la fonction setOntology sera ecrite
    return SOFA_EXON;
}

//--  --  --  --  --  --  --  --  --  --  --  --  --  --  --  --  --  --  --

/** je passe les elements un a un en affectant un vecteur de string
  *  si sofa_name est vrai (default), ce sont les noms SOFA qui sortiront, le
  *  code sinon.
  *  \todo  Les elements qui doivent �tre renseignes ne devrait pas sortir "."
  */
void
  Gff3Line::print(std::ofstream& out, bool sofa_name) const
{
    //le vecteur qui est utilis� pour l'impression
  std::vector<std::string> tab(NB_ELEMENTS_GFF3);
  tab[SEQ_ID] = (seq_id_.empty()) ? "." : seq_id_;
  tab[SOURCE] = (source_.empty()) ? "." : source_;
  if(type_sofa_==INDEFINI)
    tab[TYPE] = ".";
  else if (sofa_name)
    tab[TYPE] = Sofa::getName(type_sofa_);
  else//on fait apparaitre les zeros au debut (code sur 7 chiffres)
    tab[TYPE] = Sofa::codeToString(type_sofa_);
  //start et end peuvent-ils etre indefini ?
//  if(start_>end_)  throw exception ...
  tab[START] = (start_==INDEFINI) ? "." : to_string(start_);
  tab[END] = (end_==INDEFINI) ? "." : to_string(end_);
  tab[SCORE] = (score_==INDEFINI) ? "." : to_string(score_);
  tab[STRAND] = to_string(strand_);
  tab[PHASE] = (phase_==INDEFINI) ? "." : to_string(phase_);
  tab[ATTRIBUTES] = (attr_==0) ? "." : attr_->str();
  //ecriture
  for(short i = 0; i< NB_ELEMENTS_GFF3; ++i)
    out << tab[i] << "\t";
  out << std::endl;
}
