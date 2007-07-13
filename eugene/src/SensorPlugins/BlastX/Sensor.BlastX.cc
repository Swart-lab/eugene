// ------------------------------------------------------------------
// Copyright (C) 2004 INRA <eugene@ossau.toulouse.inra.fr>
//
// This program is open source; you can redistribute it and/or modify
// it under the terms of the Artistic License (see LICENSE file).
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
//
// You should have received a copy of Artistic License along with
// this program; if not, please see http://www.opensource.org
//
// $Id$
// ------------------------------------------------------------------
// File:     Sensor.BlastX.cc
// Contents: Sensor BlastX
// Blastx against protein databases, we simply
// enhance coding probability according to phase Dangerous
// (NR contains translation of frameshifted sequences) another
// possibility would be too forbid introns/intergenic states
// 10 levels of confidence may be used.
// ------------------------------------------------------------------

#include "Sensor.BlastX.h"

extern Parameters PAR;

/******************************************************************
 **                          SensorBlastX                        **
 *****************************************************************/

int HitsCompareLex(const void *A, const void *B)
{
    Hits **UA,**UB;

    UA = (Hits **) A;
    UB = (Hits **) B;

    if ((*UA)->Start < (*UB)->Start)
        return -1;
    if ((*UA)->Start > (*UB)->Start)
        return 1;
    if ((*UA)->NGaps > (*UB)->NGaps)
        return -1;
    if ((*UA)->NGaps < (*UB)->NGaps)
        return 1;
    if ((*UA)->Length > (*UB)->Length)
        return -1;
    if ((*UA)->Length < (*UB)->Length)
        return 1;
    return 0;
}

int HitsCompareSup(const void *A, const void *B)
{
    Hits **UA,**UB;

    UA = (Hits **) A;
    UB = (Hits **) B;

    if ((*UA)->Support > (*UB)->Support)
        return -1;
    if ((*UA)->Support < (*UB)->Support)
        return 1;
    return 0;
}

// ----------------------
// Default constructor.
// ----------------------
SensorBlastX :: SensorBlastX (int n, DNASeq *X) : Sensor(n)
{
    int i,k;
    char   tempname[FILENAME_MAX+1];
    FILE * fblast;

    N            = n;
    type         = Type_Content;
    HitTable     = NULL;
    sloppy       = PAR.getI("EuGene.sloppy");
    ppNumber     = PAR.getI("BlastX.PPNumber",N);
    stepid       = PAR.getI("Output.stepid");
    minIn        = PAR.getI("BlastX.minIn");
    levels       = PAR.getC("BlastX.levels",  N);
    intronlevels = PAR.getC("BlastX.activegaps",  N);

    Hits   *AllProt = NULL;
    NumProt = 0;

    fprintf(stderr,"Reading BlastX data, level...");
    fflush(stderr);
    inputFormat_ = to_string(PAR.getC("BlastX.format", GetNumber(),1));

    for (k = 0; k < (int)strlen(levels); k++)
    {
        // for each level given in arg:
        // e.g. with -b09 : 1st: level 0 (k=0), 2nd: level 9 (k=1)
        strcpy(tempname,PAR.getC("fstname"));
        strcat(tempname,".blast");
        i = strlen(tempname);
        tempname[i]   = levels[k];
        tempname[i+1] = 0;
	
	
	
	if ( inputFormat_ == "GFF3" )
	{
	  strcat(tempname,".gff3");
	  char * filenameSoTerms = PAR.getC("Gff3.SoTerms", GetNumber(),1);
	  char * soTerms = new char[FILENAME_MAX+1];
	  strcpy(soTerms , PAR.getC("eugene_dir"));
	  strcat(soTerms , filenameSoTerms );

	  GeneFeatureSet * geneFeatureSet = new GeneFeatureSet (tempname, soTerms);
          //geneFeatureSet->printFeature();
	  fprintf(stderr,"%c ",levels[k]);
	  fflush(stderr);
	  
	  AllProt = AllProt->ReadFromGeneFeatureSetIt(*geneFeatureSet, &NumProt, (levels[k] - '0'), 20, X);
	  
	}
	else
	{
	  // check the corresponding .blastN file (N=level given in arg)
	  fblast = FileOpen(NULL, tempname, "r", PAR.getI("EuGene.sloppy"));
  
	  if (fblast)
	  {
	      AllProt = AllProt->ReadFromFile(fblast, &NumProt, (levels[k] - '0'), 20, X->SeqLen);
	      fprintf(stderr,"%c ",levels[k]);
	      fflush(stderr);
	      fclose(fblast);
	  }
	}
	
    }
    fprintf(stderr,"done\n");
    fflush(stderr);
    Hits *ThisProt = AllProt;
    HitTable = new Hits *[NumProt+1];
    for (i = 0;  i < NumProt;  i++, ThisProt = ThisProt->Next)
    {
        HitTable[i] = ThisProt;
    }
    //Print (tempname);

    //for (i=0;  i<NumProt;  i++)
    //printf("Name:%s\tLevel:%d\n", HitTable[i]->Name,HitTable[i]->Level);
}

// ----------------------
//  Default destructor.
// ----------------------
SensorBlastX :: ~SensorBlastX ()
{
    vPos.clear();
    vPMatch.clear();
    vPMPhase.clear();
    if(HitTable != NULL)
        delete [] HitTable;
    HitTable = NULL;
}

// ----------------------
//  Init blastX.
// ----------------------
void SensorBlastX :: Init (DNASeq *X)
{
    int    i, j;
    int    Len = X->SeqLen;
    int    level, levelidx, Pphase = 0;
    float GlobalScore,    PGlobalScore = 0;

    vPos.clear();
    vPMatch.clear();
    vPMPhase.clear();

    keyBXLevel[0] = PAR.getD("BlastX.level0*", N, sloppy);
    keyBXLevel[1] = PAR.getD("BlastX.level1*", N, sloppy);
    keyBXLevel[2] = PAR.getD("BlastX.level2*", N, sloppy);
    keyBXLevel[3] = PAR.getD("BlastX.level3*", N, sloppy);
    keyBXLevel[4] = PAR.getD("BlastX.level4*", N, sloppy);
    keyBXLevel[5] = PAR.getD("BlastX.level5*", N, sloppy);
    keyBXLevel[6] = PAR.getD("BlastX.level6*", N, sloppy);
    keyBXLevel[7] = PAR.getD("BlastX.level7*", N, sloppy);
    keyBXLevel[8] = PAR.getD("BlastX.level8*", N, sloppy);
    keyBXLevel[9] = PAR.getD("BlastX.level9*", N, sloppy);
    blastxM       = PAR.getI("BlastX.blastxM*",N);
    ProtMatch      = new float[Len+1];
    ProtMatchLevel = new float[Len+1];
    ProtMatchPhase = new int[Len+1];

    for (i=0; i<= Len; i++)
    {
        ProtMatch[i]      = 0.0;
        ProtMatchLevel[i] = 0.0;
        ProtMatchPhase[i] = 0;
    }

    for (i=0;  i<NumProt;  i++)
    {

        Block *MyHSP = HitTable[i]->Match;
        level = HitTable[i]->Level;
        levelidx = rindex(levels,'0'+level)-levels;
        bool activegap = (rindex(intronlevels,'0'+level) != NULL);

        while (MyHSP != NULL)
        {
            GlobalScore = (float)(MyHSP->Score) / (float)abs(MyHSP->End-MyHSP->Start);

            // GAPS -> INTRONS ; the "intron" consideration between 2 HSP
            // requires close prot boundaries (blastxM param)
            if ( MyHSP->Prev &&
                    (abs(MyHSP->Prev->LEnd-MyHSP->LStart) <= blastxM ))
            {

                // INTRON
                if ((MyHSP->Start-MyHSP->Prev->End) >= minIn)
                {

                    for(j = MyHSP->Prev->End+1; j < MyHSP->Prev->End+1+(minIn)/2; j++)
                    {
                        // begining of the intron only, during a short region
                        // the score depends on the adjacent exon
                        if (keyBXLevel[level] >= ProtMatchLevel[j])
                        {
                            if (keyBXLevel[level] > ProtMatchLevel[j])
                            {
                                ProtMatchLevel[j]= keyBXLevel[level];
                                ProtMatch[j]= -PGlobalScore;
                                ProtMatchPhase[j] = (MyHSP->Phase > 0 ? 4 : -4);
                            }
                            else
                            {
                                if (PGlobalScore >= fabs(ProtMatch[j]))
                                {
                                    ProtMatch[j]= -PGlobalScore;
                                    ProtMatchPhase[j] = (MyHSP->Phase > 0 ? 4 : -4);
                                }
                            }
                        }
                    }

                    for (j = MyHSP->Start+1 - minIn/2; j < MyHSP->Start+1; j++)
                    {
                        // ...and the end of the intron
                        if (keyBXLevel[level] >= ProtMatchLevel[j])
                        {
                            if (keyBXLevel[level] > fabs(ProtMatchLevel[j]))
                            {
                                ProtMatchLevel[j]= keyBXLevel[level];
                                ProtMatch[j]= -GlobalScore;
                                ProtMatchPhase[j] = (MyHSP->Phase > 0 ? 4 : -4);
                            }
                            else
                                if (GlobalScore >= fabs(ProtMatch[j]))
                                {
                                    ProtMatch[j]= -GlobalScore;
                                    ProtMatchPhase[j] = (MyHSP->Phase > 0 ? 4 : -4);
                                }
                        }
                    }
                }

                // The following piece of code will take into acccount "gaps"
                // between HSP and penalise UTR/UTR introns/InterGenic
                // later. Because of protein modularity, this gives some
                // unexpected behavior and is activable per proteic bank.
                // The penalty enforced is also related to bordering
                // HSP lengths. It probably (?) would be better using the GAP
                // length.


                // INTERG: to separate InterG from Intron, 4/-4 is used to represent intron phase,
                // 0 for interG phase
                if (activegap)
                {
                    int from,to;
                    if ((MyHSP->Start-MyHSP->Prev->End) >= minIn)
                    {
                        from = MyHSP->Prev->End+1+(minIn)/2;
                        to   = MyHSP->Start - minIn/2;
                    }
                    else
                    {
                        from = MyHSP->Prev->End+1;
                        to = MyHSP->Start;
                    }

                    for (j = from; j <= to; j++)
                    {
                        // ...and the intergenic regions
                        if (keyBXLevel[level] >= ProtMatchLevel[j])
                        {
                            if (keyBXLevel[level] > fabs(ProtMatchLevel[j]))
                            {
                                ProtMatchLevel[j]= keyBXLevel[level];
                                ProtMatch[j]= -(Min(PGlobalScore,GlobalScore));
                                ProtMatchPhase[j] = 0;
                            }
                            else
                                if (GlobalScore >= fabs(ProtMatch[j]))
                                {
                                    ProtMatch[j]= -(Min(PGlobalScore,GlobalScore));
                                    ProtMatchPhase[j]= 0;
                                }
                        }
                    }
                }

                if (PAR.getI("Output.graph") && levelidx <3)
                    PlotBlastGap(MyHSP->Prev->End,Pphase,MyHSP->Start,MyHSP->Phase,levelidx);
            }

            // HITS -> CODING
            if (PAR.getI("Output.graph") && levelidx<3)
                PlotBlastHit(MyHSP->Start-1,MyHSP->End-1,MyHSP->Phase,levelidx);

            for (j = MyHSP->Start; j < MyHSP->End+1; j++)
                if (keyBXLevel[level] >= ProtMatchLevel[j])
                    if (keyBXLevel[level] > ProtMatchLevel[j])
                    {
                        ProtMatchLevel[j] = keyBXLevel[level];
                        ProtMatch[j]= GlobalScore;
                        ProtMatchPhase[j]= MyHSP->Phase;
                    }
                    else
                        if (GlobalScore >= fabs(ProtMatch[j]))
                        {
                            ProtMatch[j] = GlobalScore;
                            ProtMatchPhase[j]= MyHSP->Phase;
                        }


            PGlobalScore = GlobalScore;
            Pphase = MyHSP->Phase;
            MyHSP = MyHSP->Next;
        }
    }

    for (i = 0; i<= Len; i++)
        if(ProtMatch[i] != 0.0)
        {
            vPos.push_back    ( i );
            vPMatch.push_back ( ProtMatch[i]*ProtMatchLevel[i] );
            vPMPhase.push_back( ProtMatchPhase[i] );
        }

    delete [] ProtMatch;
    delete [] ProtMatchLevel;
    delete [] ProtMatchPhase;

    index = 0;
}


// --------------------------
//  GiveInfo Content BlastX.
// --------------------------
void SensorBlastX :: GiveInfo(DNASeq *X, int pos, DATA *d)
{
    int i;

    if ((index != 0                &&  vPos[index-1] >= pos) ||
            (index < (int)vPos.size()  &&  vPos[index]   <  pos))
    {

        iter = lower_bound(vPos.begin(), vPos.end(), pos);

        if (*iter == pos)
        {

            // EXON: phase is non 0 and different then penalize
            for (i = DATA::ExonF1; i<= DATA::ExonR3; i++)       //exons
                if (vPMPhase[iter-vPos.begin()] && (vPMPhase[iter-vPos.begin()] != PhaseAdapt(i)))
                    //	if(vPMatch[iter-vPos.begin()] < 0 ||
                    //	   ((vPMatch[iter-vPos.begin()] > 0) &&
                    //	    (vPMPhase[iter-vPos.begin()] != PhaseAdapt(i))))
                    d->contents[i] += -fabs(vPMatch[iter-vPos.begin()]);

            // always penalize these
            for (i= DATA::InterG; i<= DATA::IntronUTRR; i++)      //inter & UTRs & UTR introns
                if (vPMatch[iter-vPos.begin()] != 0)
                    d->contents[i] +=  -fabs(vPMatch[iter-vPos.begin()]);

            // penalize only vs hits
            for (i=DATA::IntronF; i<=DATA::IntronR; i++)       //introns
                if (vPMatch[iter-vPos.begin()] > 0)
                    d->contents[i] +=  -vPMatch[iter-vPos.begin()];
            index = iter-vPos.begin() + 1;

        }
        else
            index = iter-vPos.begin();
    }
    else if( index < (int)vPos.size()  &&  vPos[index] == pos )
    {
        // EXON: phase is non 0 and different then penalize
        for (i = DATA::ExonF1; i<= DATA::ExonR3; i++)       //exons
            if (vPMPhase[index] && (vPMPhase[index] != PhaseAdapt(i)))
                //if(vPMatch[index] < 0 ||
                //	 ((vPMatch[index] > 0) && (vPMPhase[index] != PhaseAdapt(i))))
                d->contents[i] += -fabs(vPMatch[index]);

        // always penalize these
        for (i= DATA::InterG; i<= DATA::IntronUTRR; i++)      //inter & UTRs & UTR introns  !! mettre IntronUTRR
            if (vPMatch[index] != 0)
                d->contents[i] +=  -fabs(vPMatch[index]);

        // penalize only vs hits
        for (i=DATA::IntronF; i<=DATA::IntronR; i++)       //introns
            if (vPMatch[index] > 0)
                d->contents[i] +=  -vPMatch[index];
        index++;
    }
}

// ----------------------------
//  Plot Sensor information
// ----------------------------
void SensorBlastX :: Plot(DNASeq *X)
{
    //  Plot done in init
}

// ------------------
//  Post analyse
// ------------------
void SensorBlastX :: PostAnalyse(Prediction *pred, FILE *MINFO)
{
    int state    = 0, start  = 0, end = 0;
    int cdsStart = 0, cdsEnd = 0;
    int CodingNuc    = 0;
    int SupportedNuc = 0;

    int pprocess = PAR.getI("BlastX.PostProcess",GetNumber());
    if (!pprocess)
        return;
    if (NumProt == 0)
        return;

    fprintf(MINFO, "#=========================================#\n");
    fprintf(MINFO, "#=           Protein evidences           =#\n");
    fprintf(MINFO, "#=========================================#\n");

    index=0;
    qsort((void *)HitTable, NumProt, sizeof(void *), HitsCompareLex);

    // Reset static in Prot Support or FEA Support
    ProtSupport(NULL, NULL, 100, 0, NULL, 0, 0);

    for (int i=0; i<pred->nbGene; i++)
    {
        SupportedNuc = 0;
        CodingNuc    = pred->vGene[i]->exLength;
        cdsStart     = pred->vGene[i]->cdsStart + 1;
        cdsEnd       = pred->vGene[i]->cdsEnd   + 1;

        if (pprocess == 1)
        {

            // WARNING: this analysis does only take into account the hits
            // that have not be "shadowed" by other stronger hits in the
            // Init phase.

            for (int j=0; j<pred->vGene[i]->nbFea(); j++)
            {
                state = pred->vGene[i]->vFea[j]->state;
                start = pred->vGene[i]->vFea[j]->start;
                end   = pred->vGene[i]->vFea[j]->end;
                if (State2Status[state] == 3)
                { //coding
                    // index=indice du premier match etant > ou = au debut de l'exon
                    while( (index < (int)vPos.size()) && (vPos[index] < start-1) )
                        index++;

                    if (index < (int)vPos.size())
                    { // il reste des hits
                        // pour chaque nuc de l'exon supporte par un hit
                        while (index < (int)vPos.size() && vPos[index]<end)
                        {
                            if (State2Phase[state] ==  vPMPhase[index])
                                SupportedNuc++;
                            index++;
                        }
                    }
                }
            }
            fprintf(MINFO, "      CDS (prot.)  %7d %7d    %5d     ",
                    cdsStart, cdsEnd, cdsEnd - cdsStart + 1);
            fprintf(MINFO, "supported on %d/%d coding.\n", SupportedNuc, CodingNuc);
        }
        else
            ProtSupport(pred, MINFO, cdsStart, cdsEnd, HitTable, NumProt, i+1);
    }
}

// -------------------------------------------------------------------------
//  Post analyse (Analyse de la prédiction (features) par rapport aux Prot)
//    debut/fin  = debut/fin de traduit
// -------------------------------------------------------------------------
void SensorBlastX :: ProtSupport(Prediction *pred, FILE *MINFO, int debut,
                                 int fin, Hits **HitTable, int Size, int NumGene)
{
    static int ProtIndex;
    Block *ThisBlock;
    int SupportProt, i, j;
    int from  = 0, to = 0;
    int state = 0;
    std::vector <int> vSupProtI;                // index prot supportant pred
    Hits **TMPHitTable = new Hits *[NumProt+1]; // need for sort by % support

    if (pred == NULL)
    {
        ProtIndex = 0;
        return;
    }

    /***********************************************************************/
    /** Objectif : obtenir un vecteur contenant les index des prots       **/
    /**            qui supportent la prédiction                           **/
    /***********************************************************************/
    // si l'iteration precedente a atteint l'extremite
    if (ProtIndex >= Size)
        ProtIndex = Max(0, Size-1);

    // on rembobine....
    while ((ProtIndex > 0) && (HitTable[ProtIndex]->End > debut))
        ProtIndex--;

    if (ProtIndex >= 0  &&  HitTable[ProtIndex]->End < debut)
        ProtIndex++;

    while (ProtIndex >=0  &&  ProtIndex < Size)
    {
        // la derniere prot exploitable est passee
        if (HitTable[ProtIndex]->Start > fin)
            break;

        SupportProt = 0;
        ThisBlock = HitTable[ProtIndex]->Match;

        while (ThisBlock)
        {
            from = Max(debut, ThisBlock->Start);
            to   = Min(fin,   ThisBlock->End);
            if (from < to)
            {
                SupportProt = 1;
                break;
            }
            ThisBlock = ThisBlock->Next;
        }

        // Si Prot supporte potentiellement
        if (SupportProt)
            vSupProtI.push_back( ProtIndex );
        ProtIndex++;
    }


    /***********************************************************************/
    /* Objectif : Analyser chaque feature predite -> supportée ?           */
    /***********************************************************************/
    int  start = 0, end = 0, len;
    char fea[5];
    char strand = pred->vGene[NumGene-1]->vFea[0]->strand;
    int  codingNuc = 0;

    for (i=0; i<pred->vGene[NumGene-1]->nbFea(); i++)
    {
        state = pred->vGene[NumGene-1]->vFea[i]->state;
        start = pred->vGene[NumGene-1]->vFea[i]->start;
        end   = pred->vGene[NumGene-1]->vFea[i]->end;

        if ( (start >= debut)  &&  (end <= fin))
        {
            len      = 0;
            int numF = -1;

            if (state <= TermR3  ||  (state >= UTR5F && state <= UTR3R))
                numF = pred->vGene[NumGene-1]->vFea[i]->number;


            if (state <= TermR3)
            {
                strcpy(fea, "Exon");
 				codingNuc += (end - start + 1);
            }
            else if(state == UTR5F  ||  state == UTR5R)
                strcpy(fea, "UTR5");
            else if(state == UTR3F  ||  state == UTR3R)
                strcpy(fea, "UTR3");

            if (numF != -1)
            {
                if ((int)vSupProtI.size() > 0)
                    // Longueur totale supportée par les prots
                    len = LenSup(HitTable, pred, vSupProtI, -1, start, end);

                if (len > 0)
                {
                    fprintf(MINFO, "%s.%d.%d\tEuGene_prot\t%s\t%d\t%d\t%d\t%c\t.\t",
                            pred->seqName, (((NumGene-1)*stepid)+1), numF,
                            fea,start,end,len,strand);
                    fprintf(MINFO, "%d\t", (int)((float)len/(end-start+1)*100));

                    for (int k=0;  k<NumProt;  k++)
                        HitTable[k]->Support = 0;

                    for (j=0; j<(int)vSupProtI.size(); j++)
                    {
                        // Longueur supportée par la prot
                        len = LenSup(HitTable, pred, vSupProtI, j, start, end);
                        HitTable[vSupProtI[j]]->Support = (int)((float)len/(end-start+1)*100);
                    }
                    // On copie la hittable pour trier sur le % supporté
                    for (int k=0; k<NumProt; k++)
                        TMPHitTable[k] = HitTable[k];
                    qsort((void*)TMPHitTable, NumProt, sizeof(void*), HitsCompareSup);

                    // On affiche les ppNumber premiers hits supportant
                    for (j=0; j<NumProt && j<ppNumber && TMPHitTable[j]->Support!=0;j++)
                        fprintf(MINFO, "%s(%d,%d) ", TMPHitTable[j]->Name,
                                TMPHitTable[j]->Support, TMPHitTable[j]->Level);
                    fprintf(MINFO, "\n");
                }
            }
        }
    }

    /***********************************************************************/
    /* Objectif : Analyser la CDS et le gene predit -> supportée ?         */
    /***********************************************************************/
    start = debut;
    end   = fin;
    strcpy(fea, "CDS");

    if (end >= start)
    {
        len = 0;
        if((int)vSupProtI.size() > 0)
            // Longueur totale supportée par les prots
            len = LenSup(HitTable, pred, vSupProtI, -1, start, end);

        if(len > 0)
        {
            fprintf(MINFO, "%s.%d  \tEuGene_prot\t%s\t%d\t%d\t%d\t%c\t.\t",
                    pred->seqName, (((NumGene-1)*stepid)+1), fea,
                    start, end, len, strand);
            fprintf(MINFO, "%d\t", (int)((float)len/codingNuc*100));
            for (int k=0;  k<NumProt;  k++)
                HitTable[k]->Support = 0;

			if ( (int)((float)len/codingNuc*100) > 100 )
			{
				exit(3);
			}

            for(j=0; j<(int)vSupProtI.size(); j++)
            {
                // Longueur supportée par la prot
                len = LenSup(HitTable, pred, vSupProtI, j, start, end);
                HitTable[vSupProtI[j]]->Support = (int)((float)len/codingNuc*100);
            }
            // On copie la hittable pour trier sur le % supporté
            for (int k=0; k<NumProt; k++)
                TMPHitTable[k] = HitTable[k];
            qsort((void*)TMPHitTable, NumProt, sizeof(void*), HitsCompareSup);

            // On affiche les ppNumber premiers hits supportant
            for(j=0; j<NumProt && j<ppNumber && TMPHitTable[j]->Support!=0; j++)
                fprintf(MINFO, "%s(%d,%d) ", TMPHitTable[j]->Name,
                        TMPHitTable[j]->Support, TMPHitTable[j]->Level);
            fprintf(MINFO, "\n");
        }
    }
    vSupProtI.clear();
    if(TMPHitTable != NULL)
        delete [] TMPHitTable;
    TMPHitTable = NULL;
    return;
}

// -------------------------------------------------------------------------
//  Length supported by Prot.
//    If index = -1 -> Length supported by ALL Prot.
//    Else          -> Length supported by ONE Prot.
// -------------------------------------------------------------------------
int SensorBlastX :: LenSup(Hits **HitTable, Prediction *pred,
                           std::vector<int> vSupProtI,
                           int index, int beg, int end)
{
    int supported = 0;
//    unsigned char *Sup;
	short *Sup;
    Block *ThisBlock;
    int from  = 0, to = 0;
    int state = 0;
    int i, j;

    // Sup = new unsigned char[end-beg+1];
	Sup = new short[end-beg+1];
    for (i=0; i<=(end-beg); i++)
        Sup[i]=0;

    i = (index == -1  ?  0 : index);

    for(; i<(int)vSupProtI.size() || i==index; i++)
    {
        ThisBlock = HitTable[vSupProtI[i]]->Match;
        while (ThisBlock)
        {
            from = Max(beg, ThisBlock->Start+1);
            to   = Min(end, ThisBlock->End+1);

            for (j = from; j <= to; j++)
            {
                state = pred->GetStateForPos(j);
                if ((Sup[j-beg] ==0) &&  (State2Phase[state] == ThisBlock->Phase) && state <= TermR3 )
                {
                    Sup[j-beg] = 1;
                    supported++;
                }
            }
            ThisBlock = ThisBlock->Next;
        }
        if (index != -1)
            break;
    }
    delete [] Sup;
    return supported;
}


void SensorBlastX :: Print (char name[FILENAME_MAX+1])
{
  FILE *fp;
  strcat (name, ".out");
  if (!(fp = fopen(name, "w"))) {
    fprintf(stderr, "cannot write in %s\n",  name);
    exit(2);
  }
 for (int i=0;  i<NumProt;  i++)
  fprintf(fp, "Name:%s\tLevel:%d\n", HitTable[i]->Name,HitTable[i]->Level);
  


  fclose(fp);
}
