/* ---------------------------------------------------------------------------
* Fichier : lex.l
* Module : lex
* Version : 1.0
* Date de creation de la version : 01/07/99
* Auteurs : Lucien Duret
* ---------------------
* Description : contient l'analyseur lexicale utilisateur 
* pour FrameD++ et EuGene
* ---------------------
* Divers : 
* ---------------------------------------------------------------------------
* --- INRA --------------  Station de Biometrie et d'Intelligence Artificielle
* ---------------------------------------------------------------------------
*/

/*    ANALYSEUR LEXICAL                                          */
/*     repere les elements de la grammaire :                     */
/*       - RAISON START STOP DONNEUR ACCEPTEUR                   */
/*       - EXON INTRON INTERGENIQUE NUM CERTITUDE                     */


%{
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
// #include "yacc.tab.h" 

extern void yyerror(char* s);   /* nouveauté Linux compatible Solaris */


%}

WS [ \t]+
NL [\n]+  
RC [\r]+

%%

{WS} 
{NL} {return(NL);}
{RC}   


"#".*    {yylval.chaine = strdup(yytext);return(RAISON);}

"start"  {yylval.tab = 13;return(TYPE);}

"stop" {yylval.tab = 15; return(TYPE);}

"donor"  {yylval.tab = 19;return(TYPE);}

"acceptor"  {yylval.tab = 17;return(TYPE);}

"exon"  {yylval.tab = 0;return(EXON);}

"intron"  {yylval.tab = 0;return(INTRON);}

"utr3" {yylval.tab = 0;return(UTR3);}

"utr5" {yylval.tab = 0;return(UTR5);}

"intergenic"  {yylval.tab = 0;return(INTERGENIQUE);}

"f"  { yylval.entier = 0 ; return(BRIN);}
"r"  { yylval.entier = 1; return(BRIN);}

-?[0-9]+  {yylval.entier = atoi( yytext ); return(NUM);}

"[" {return(CO);}
"]" {return(CF);}
"."+ {return(PT);}

"check"  { yylval.check = check ; return(CH);}
"nocheck"  { yylval.check = nocheck; return(CH);}

"infinity" {yylval.reel = -log(0.0); return(POIDS);}

"-infinity" {yylval.reel = log(0.0); return(POIDS);}


%-?[0-9]|[0-9]+"."[0-9]+|-?0x[0-9a-h]"."?[0-9a-h]*p[+-][0-9]+  {yylval.reel = strtod(yytext,NULL); return(POIDS);}


-?[0-9]+"."[0-9]+  {yylval.reel = strtod(yytext,NULL); return(POIDS);}

.  printf( "Unknown character : %s\n", yytext );


%%

