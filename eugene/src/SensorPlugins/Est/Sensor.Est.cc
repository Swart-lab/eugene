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
// File:     Sensor.Est.cc
// Contents: Sensor
// ------------------------------------------------------------------

#include "Sensor.Est.h"
#include "../../MSensor.h"

#define EXTREMEMARGIN 1
#define Inconsistent(x) (((x) & Hit) && ((x) & Gap))

extern MasterSensor *MS;
extern Parameters   PAR;

/*************************************************************
 **                        SensorEst                        **
 *************************************************************/
int HitsCompare(const void *A, const void *B)
{
    Hits **UA,**UB;

    UA = (Hits **) A;
    UB = (Hits **) B;

    if ((*UA)->NGaps > (*UB)->NGaps)
        return -1;
    if ((*UA)->NGaps < (*UB)->NGaps)
        return  1;
    if ((*UA)->Length > (*UB)->Length)
        return -1;
    if ((*UA)->Length < (*UB)->Length)
        return  1;
    return strcmp((*UA)->Name,(*UB)->Name); //pour un tri sans ambiguité
}

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
//  Default constructor.
// ----------------------
SensorEst :: SensorEst (int n, DNASeq *X) : Sensor(n)
{
    FILE *fEST;
    char tempname[FILENAME_MAX+1];
    int  i;

    type     = Type_Content;
    HitTable = NULL;
    N        = n;

    vPos.clear();
    vESTMatch.clear();

    estM           = PAR.getI("Est.estM",N);
    utrM           = PAR.getI("Est.utrM",N);
    ppNumber       = PAR.getI("Est.PPNumber",N);
    stepid         = PAR.getI("Output.stepid",N);
    MinDangling    = PAR.getI("Est.MinDangling",N);
    MaxIntron      = PAR.getI("Est.MaxIntron",N);
    MaxIntIntron   = PAR.getI("Est.MaxInternalIntron",N);
    DonorThreshold = PAR.getD("Est.StrongDonor",N);
    DonorThreshold = log(DonorThreshold/(1-DonorThreshold));

    index = 0;

    ESTMatch = new unsigned char[X->SeqLen+1];
    for (i = 0; i <= X->SeqLen; i++)
        ESTMatch[i] = 0;

    fprintf(stderr,"Reading cDNA hits............");
    fflush(stderr);

    strcpy(tempname, PAR.getC("fstname"));
    strcat(tempname, ".est");
    NumEST = 0;
    Hits * AllEST=NULL;
      
    inputFormat_ = to_string(PAR.getC("Est.format", GetNumber(),1));

    if ( inputFormat_ == "GFF3" )
    {
      strcat(tempname,".gff3");
      char * filenameSoTerms = PAR.getC("Gff3.SoTerms", GetNumber(),1);
      char * soTerms = new char[FILENAME_MAX+1];
      strcpy(soTerms , PAR.getC("eugene_dir"));
      strcat(soTerms , filenameSoTerms );

      GeneFeatureSet * geneFeatureSet = new GeneFeatureSet (tempname, soTerms);
      geneFeatureSet->printFeature();
      AllEST = AllEST->ReadFromGeneFeatureSetIt(*geneFeatureSet, &NumEST, -1, 0, X);
      HitTable = ESTAnalyzer(AllEST, ESTMatch, estM, &NumEST, X);
    }
    else
    {
      fEST = FileOpen(NULL, tempname, "r", PAR.getI("EuGene.sloppy"));
      if (fEST)
      {   //ReadFromFile (EstFile  EstNumber  Level  Margin)
	  AllEST = AllEST->ReadFromFile(fEST, &NumEST, -1, 0,X->SeqLen);
	  HitTable = ESTAnalyzer(AllEST, ESTMatch, estM, &NumEST, X);
	  fclose(fEST);
      }
      else
      {
	  fprintf(stderr,"Error while openning %s\n",tempname);
	  fflush(stderr);
      }
    }
    
    for (i = 0; i<= X->SeqLen; i++)
        if(ESTMatch[i] != 0)
        {
            vPos.push_back      ( i );
            vESTMatch.push_back ( ESTMatch[i] );
        }
    //Print(tempname);
    delete [] ESTMatch;
}

// ----------------------
//  Default destructor.
// ----------------------
SensorEst :: ~SensorEst ()
{
    vPos.clear();
    vESTMatch.clear();

    if(HitTable != NULL)
        delete [] HitTable;
    HitTable = NULL;
}

// --------------------------------
//  Init est.
//  Exploiting spliced alignements
//  against EST and complete cDNA.
// --------------------------------
void SensorEst :: Init (DNASeq *X)
{
    estP = PAR.getD("Est.estP*",N);
    utrP = PAR.getD("Est.utrP*",N);
    spliceBoost    = PAR.getD("Est.SpliceBoost*",N);
    //for(int jj=0;jj<(int)vPos.size();jj++)
    //printf("vPos[%d]:%d\tvESTM[%d]:%d\n",jj,vPos[jj]+1,jj,vESTMatch[jj]);
}

// -----------------------
//  GiveInfo signal est.
// -----------------------
void SensorEst :: GiveInfo (DNASeq *X, int pos, DATA *d)
{
    unsigned char cESTMatch = 0; // current ESTMatch

    // Peut on faire un bete acces sequentiel ?
    if((index != 0                &&  vPos[index-1] >= pos) ||
            (index < (int)vPos.size()  &&  vPos[index]   <  pos))
    {
        // Non... on se repositionne en temps logarithmique
        iter = lower_bound(vPos.begin(), vPos.end(), pos);
        index = iter-vPos.begin();
    }
    // On est juste avant ou sur pos

    // Si on est dessus
    if (index < (int)vPos.size()  &&  vPos[index] == pos)
    {

        cESTMatch = vESTMatch[index];

        // Favor splice sites in marginal exon-intron regions
        if ((cESTMatch & MLeftForward) &&
                d->sig[DATA::Don].weight[Signal::Forward] != 0.0)
            d->sig[DATA::Don].weight[Signal::Forward] += spliceBoost;

        if ((cESTMatch & MRightForward) &&
                d->sig[DATA::Acc].weight[Signal::Forward] != 0.0)
            d->sig[DATA::Acc].weight[Signal::Forward] += spliceBoost;

        if ((cESTMatch & MLeftReverse) &&
                d->sig[DATA::Acc].weight[Signal::Reverse])
            d->sig[DATA::Acc].weight[Signal::Reverse] += spliceBoost;

        if ((cESTMatch & MRightReverse) &&
                d->sig[DATA::Don].weight[Signal::Reverse] != 0.0)
            d->sig[DATA::Don].weight[Signal::Reverse] += spliceBoost;

        // Exon ou UTR  Forward
        // Si on a un Gap EST ou si l'on connait le sens du match EST
        if ((cESTMatch & Gap) ||
                ((cESTMatch & Hit) && !(cESTMatch & HitForward)))
        {
            d->contents[DATA::UTR5F] += estP;
            d->contents[DATA::UTR3F] += estP;
            for(int i=0; i<3; i++)
                d->contents[i] += estP;
        }

        // Exon ou UTR Reverse
        // Si on a un Gap EST ou si l'on connait le sens du match EST
        if ((cESTMatch & Gap) ||
                ((cESTMatch & Hit) && !(cESTMatch & HitReverse)))
        {
            d->contents[DATA::UTR5R] += estP;
            d->contents[DATA::UTR3R] += estP;
            for(int i=3; i<6; i++)
                d->contents[i] += estP;
        }

        // Intron Forward
        // Si on a un Hit EST ou si l'on connait le sens du match EST
        if((cESTMatch & Hit) ||
                ((cESTMatch & Gap) && !(cESTMatch & GapForward)))
        {
            d->contents[DATA::IntronF] += estP;
            d->contents[DATA::IntronUTRF] += estP;
        }

        // Intron Reverse
        // Si on a un Hit EST ou si l'on connait le sens du match EST
        if((cESTMatch & Hit) ||
                ((cESTMatch & Gap) && !(cESTMatch & GapReverse)))
        {
            d->contents[DATA::IntronR] += estP;
            d->contents[DATA::IntronUTRR] += estP;
        }

        // Intergenique: tout le temps si on a un match
        d->contents[DATA::InterG] += ((cESTMatch & (Gap|Hit)) != 0)*estP;

        d->EstMatch = cESTMatch;  // WARNING : EST -> on est dans intron
        index++;
    }

    // Pour que les UTR soient supportés par un EST
    if (cESTMatch == 0  &&  (int)vPos.size() != 0)
    {
        // Left
        for (int k=1; k<=utrM; k++)
        {
            iter = lower_bound(vPos.begin(), vPos.end(), pos+k);
            if(*iter == pos+k)
            {
                cESTMatch = vESTMatch[iter-vPos.begin()];
                // If only Margin (-> EST extremities) then penalize all utr tracks
                if((cESTMatch & Margin) && !(cESTMatch & Gap) && !(cESTMatch & Hit))
                {
                    d->contents[DATA::UTR5F] += log(utrP);
                    d->contents[DATA::UTR5R] += log(utrP);
                    d->contents[DATA::UTR3F] += log(utrP);
                    d->contents[DATA::UTR3R] += log(utrP);
                    break;
                }
            }
        }
        // Right
        for (int k=utrM; k>0; k--)
        {
            iter = lower_bound(vPos.begin(), vPos.end(), pos-k);
            if(*iter == pos-k)
            {
                cESTMatch = vESTMatch[iter-vPos.begin()];
                // If only Margin (-> EST extremities) then penalize all utr tracks
                if((cESTMatch & Margin) && !(cESTMatch & Gap) && !(cESTMatch & Hit))
                {
                    d->contents[DATA::UTR5F] += log(utrP);
                    d->contents[DATA::UTR5R] += log(utrP);
                    d->contents[DATA::UTR3F] += log(utrP);
                    d->contents[DATA::UTR3R] += log(utrP);
                    break;
                }
            }
        }
    }
}

// -----------------------
//  ESTAnalyzer.
// -----------------------
Hits** SensorEst :: ESTAnalyzer(Hits *AllEST, unsigned char *ESTMatch,
                                int EstM, int *NumEST, DNASeq *X)
{
    int i,j,k;
    int Rejected = 0;
    Hits *ThisEST = NULL;
    Block *ThisBlock = NULL;
    
    
    fprintf(stderr,"%d sequences read\n",*NumEST);
    fflush(stderr);

    // on trie les hits sur le nombre de gaps et la
    // longueur. L'idee est d'eliminer les epissages partiels et de
    // favoriser la longueur de toute facon.
    Hits **HitTable = new Hits *[*NumEST+1];
    for (i = 0, ThisEST = AllEST; i < *NumEST; i++, ThisEST = ThisEST->Next)
    {        HitTable[i] = ThisEST;    }
    // pour memoriser le premier (a liberer)
    HitTable[*NumEST] = AllEST;

    qsort((void *)HitTable, *NumEST, sizeof(void *), HitsCompare);

    for (int index = 0; index < *NumEST; index++)
    {
        int Inc;
        int ExonInc,TheStrand;
        double WorstSpliceF, WorstSpliceR;
        double DonF,AccF,DonR,AccR;
        DATA dTmp;

        ThisEST = HitTable[index];

        // Le veritable  brin est a priori indetermine
        TheStrand = HitForward | HitReverse;
        Inc = 0;
        ExonInc = 0;
        WorstSpliceF = WorstSpliceR = -NINFINITY;

        // Look for each match in the current Hit
        ThisBlock = ThisEST->Match;

        // First hit: delete short dangling
        if ((ThisBlock->Next != NULL) &&
                ((abs(ThisBlock->Start - ThisBlock->End) <= MinDangling) ||
                 (abs(ThisBlock->End - ThisBlock->Next->Start) >= MaxIntron)))
        {
            ThisEST->Match = ThisBlock->Next;
            ThisBlock->Next->Prev = NULL;
            fprintf(stderr,
                    "   [%s]: Suspicious dangling match removed (len. %d, gap %d)\n",
                    ThisEST->Name,abs(ThisBlock->Start - ThisBlock->End),
                    abs(ThisBlock->End - ThisBlock->Next->Start));
            ThisBlock->Next = NULL;
            delete ThisBlock;
            ThisEST->NGaps--;
            ThisBlock = ThisEST->Match;
        }

        // First Step: tries to determine strandedness
        while (ThisBlock)
        {

            // Last block: delete short dangling
            if ((ThisBlock->Next == NULL) && (ThisBlock->Prev != NULL) &&
                    ((abs(ThisBlock->Start-ThisBlock->End) <= MinDangling) ||
                     (abs(ThisBlock->Start-ThisBlock->Prev->End) >= MaxIntron)))
            {
                ThisBlock->Prev->Next = NULL;
                fprintf(stderr,
                        "   [%s]: Suspicious dangling match removed (len. %d, gap %d)\n",
                        ThisEST->Name,abs(ThisBlock->Start-ThisBlock->End),
                        abs(ThisBlock->Start-ThisBlock->Prev->End));
                delete  ThisBlock;
                ThisBlock = NULL;
                ThisEST->NGaps--;
                break;
            }

            // si on a un gap ?
            if ((ThisBlock->Prev != NULL) &&
                    abs(ThisBlock->Prev->LEnd - ThisBlock->LStart) <= 6)
            {
                DonF = NINFINITY;
                DonR = NINFINITY;
                AccF = NINFINITY;
                AccR = NINFINITY;

                //	printf("[%d %d]\n",ThisBlock->Start,ThisBlock->End);

                for (j = -EstM; j <= EstM; j++)
                {
                    k = Min(X->SeqLen,Max(0,ThisBlock->Prev->End+j+1));
                    MS->GetInfoSpAt(Type_Acc|Type_Don, X, k, &dTmp);

                    DonF = Max(DonF, dTmp.sig[DATA::Don].weight[Signal::Forward]-
                               dTmp.sig[DATA::Don].weight[Signal::ForwardNo]);
                    AccR = Max(AccR, dTmp.sig[DATA::Acc].weight[Signal::Reverse]-
                               dTmp.sig[DATA::Acc].weight[Signal::ReverseNo]);

                    k = Min(X->SeqLen,Max(0,ThisBlock->Start+j));
                    if(MS->GetInfoSpAt(Type_Acc|Type_Don, X, k, &dTmp))
                    {
                        DonR = Max(DonR, dTmp.sig[DATA::Don].weight[Signal::Reverse]-
                                   dTmp.sig[DATA::Don].weight[Signal::ReverseNo]);
                        AccF = Max(AccF, dTmp.sig[DATA::Acc].weight[Signal::Forward]-
                                   dTmp.sig[DATA::Acc].weight[Signal::ForwardNo]);
                    }
                    else
                    {
                        fprintf(stderr,"   WARNING: cDNA hits ignored."
                                " No splices sites predicted !\n");
                        exit(2);
                    }
                }

                //	printf("Extreme splices: %f %f - %f %f\n",DonF,AccF,DonR,AccR);
                WorstSpliceF = Min(WorstSpliceF,DonF);
                WorstSpliceF = Min(WorstSpliceF,AccF);
                WorstSpliceR = Min(WorstSpliceR,DonR);
                WorstSpliceR = Min(WorstSpliceR,AccR);
            }
            ThisBlock = ThisBlock->Next;
        }

        //    printf("Extreme splices: %e %e\n",WorstSpliceF,WorstSpliceR);

        // Tous les blocs ont ete traites
        if (WorstSpliceF == NINFINITY)
            TheStrand &= (~HitForward);
        if (WorstSpliceR == NINFINITY)
            TheStrand &= (~HitReverse);

        // next iteration on the same EST
        ThisBlock = ThisEST->Match;
        while (TheStrand && ThisBlock)
        {
            // Check for consistency with already read Hits
            // The Inc flag will keep the Inconsistency status
            // 1 - inconsistent with a previous EST
            // 2 - no splice site on the borders of a gap
            // 3 - an exon contains STRONG donor on both strands

            for (i = ThisBlock->Start+EstM; !Inc && i <= ThisBlock->End-EstM; i++)
            {

                if (((ESTMatch[i] & Hit) && !(ESTMatch[i] & TheStrand)) ||
                        (Inconsistent(ESTMatch[i] | TheStrand)))
                {
                    fprintf(stderr,"   [%s]: inconsistent hit [%d-%d]\n",
                            ThisEST->Name,ThisBlock->Start+1,ThisBlock->End+1);
                    Inc = 1;
                }
            }

            // si on a un gap
            if ((ThisBlock->Prev != NULL) &&
                    abs(ThisBlock->Prev->LEnd - ThisBlock->LStart) <= 6)
            {
                for (i=ThisBlock->Prev->End+1+EstM; !Inc && i<ThisBlock->Start-EstM; i++)
                    if (((ESTMatch[i] & Gap) && !(ESTMatch[i] & (TheStrand << HitToGap))) ||
                            (Inconsistent(ESTMatch[i] | (TheStrand << HitToGap))))
                    {
                        fprintf(stderr,"   [%s]: inconsistent gap [%d-%d]\n",
                                ThisEST->Name,ThisBlock->Prev->End+2,ThisBlock->Start);
                        Inc = 1;
                    }
            }

            DonF = NINFINITY;
            DonR = NINFINITY;

            // calcul des sites d'epissage internes a l'exon
            for (i = ThisBlock->Start+EstM+1; !Inc && i <= ThisBlock->End-EstM-1; i++)
            {
                MS->GetInfoSpAt(Type_Acc|Type_Don, X, i, &dTmp);
                DonF = Max(DonF, dTmp.sig[DATA::Don].weight[Signal::Forward]-
                           dTmp.sig[DATA::Don].weight[Signal::ForwardNo]);
                DonR = Max(DonR, dTmp.sig[DATA::Don].weight[Signal::Reverse]-
                           dTmp.sig[DATA::Don].weight[Signal::ReverseNo]);
            }
            if (DonF > DonorThreshold)
                ExonInc |= 1;
            if (DonR > DonorThreshold)
                ExonInc |= 2;
            if (ExonInc == 3 && !Inc && !ThisEST->NGaps)
            {
                fprintf(stderr,"   [%s]: Gapless EST with strong donor [%d-%d]\n",
                        ThisEST->Name,ThisBlock->Start+1,ThisBlock->End+1);
                Inc = 3;
            }

            ThisBlock = ThisBlock->Next;
        }

        if (!TheStrand)
        {
            fprintf(stderr, "   [%s]: no matching splice site\n",ThisEST->Name);
            Inc =2;
        }

        // Si une incoherence est detectee, on va jeter la sequence

        if (Inc)
        {
            Rejected++;
            ThisEST->Rejected = 1;

            ThisBlock = ThisEST->Match;
            while (ThisBlock)
            {

                if (TheStrand & HitForward)
                    PlotESTHit(ThisBlock->Start,ThisBlock->End,1,1);
                if (TheStrand & HitReverse)
                    PlotESTHit(ThisBlock->Start,ThisBlock->End,-1,1);

                if (PAR.getI("Output.graph") && (ThisBlock->Prev != NULL) &&
                        abs(ThisBlock->Prev->LEnd-ThisBlock->LStart) <= 6)
                {
                    if (TheStrand & HitForward)
                        PlotESTGap(ThisBlock->Prev->End,ThisBlock->Start,1,1);
                    if (TheStrand & HitReverse)
                        PlotESTGap(ThisBlock->Prev->End,ThisBlock->Start,-1,1);
                }
                ThisBlock = ThisBlock->Next;
            }
        }
        // sinon on l'exploite
        else
        {
            int LBoundary;
            int RBoundary;
            unsigned char Info;

            ThisBlock = ThisEST->Match;

            while (ThisBlock)
            {

                // Aligners tend to extend beyond the true hit on
                // extremities: we remove EstM on frontiers

                LBoundary = ((ThisBlock->Prev == NULL) ?
                             Min(X->SeqLen, ThisBlock->Start+EstM) :
                             ThisBlock->Start);

                RBoundary = ((ThisBlock->Next == NULL) ?
                             Max(0,ThisBlock->End-EstM) :
                             ThisBlock->End);

                for (i = LBoundary; i <= RBoundary; i++)
                {
                    if ((Info = (ESTMatch[i] & Hit)))
                        ESTMatch[i] |= (Info & TheStrand);
                    else
                        ESTMatch[i] |= TheStrand;
                }

                // Do we want to mark extreme region or alignements as margins ?
#ifdef EXTREMEMARGIN
                if (ThisBlock->Prev == NULL)
                    for (i = ThisBlock->Start; i<LBoundary; i++)
                        ESTMatch[i] |= TheStrand << HitToMLeft;

                if (ThisBlock->Next == NULL)
                    for (i = RBoundary+1; i <= ThisBlock->End; i++)
                        ESTMatch[i] |= TheStrand << HitToMRight;
#endif

                if (PAR.getI("Output.graph"))
                {
                    if (TheStrand & HitForward)
                        PlotESTHit(ThisBlock->Start,ThisBlock->End,1,0);
                    if (TheStrand & HitReverse)
                        PlotESTHit(ThisBlock->Start,ThisBlock->End,-1,0);
                }

                if ((ThisBlock->Prev != NULL) &&
                        abs(ThisBlock->Prev->LEnd-ThisBlock->LStart) <= 6)
                {

                    for (i = Max(0,ThisBlock->Prev->End+1-EstM);
                            i < Min(X->SeqLen, ThisBlock->Prev->End+1+EstM); i++)
                        if ((Info = (ESTMatch[i] & MLeft)))
                            ESTMatch[i] |= (Info & (TheStrand << HitToMLeft));
                        else
                            ESTMatch[i] |= TheStrand << HitToMLeft;

                    if (abs(ThisBlock->Prev->End - ThisBlock->Start) < MaxIntIntron)
                    {

                        for (i = ThisBlock->Prev->End+1; i <= ThisBlock->Start-1; i++)
                            if ((Info = (ESTMatch[i] & Gap)))
                                ESTMatch[i] |= (Info & (TheStrand << HitToGap));
                            else
                                ESTMatch[i] |= TheStrand << HitToGap;
                    }
                    else
                        fprintf(stderr, "   [%s]: long intron ignored (%d bases)\n",
                                ThisEST->Name,abs(ThisBlock->Prev->End - ThisBlock->Start));

                    for (i =  Max(0,ThisBlock->Start-1-EstM);
                            i <  Min(X->SeqLen, ThisBlock->Start-1+EstM); i++)
                        if ((Info = (ESTMatch[i] & MRight)))
                            ESTMatch[i] |= (Info & (TheStrand << HitToMRight));
                        else
                            ESTMatch[i] |= TheStrand << HitToMRight;


                    if (PAR.getI("Output.graph"))
                    {
                        if (TheStrand & HitForward)
                            PlotESTGap(ThisBlock->Prev->End,ThisBlock->Start,1,0);
                        if (TheStrand & HitReverse)
                            PlotESTGap(ThisBlock->Prev->End,ThisBlock->Start,-1,0);
                    }
                }
                ThisBlock = ThisBlock->Next;
            }
        }
    }

    if (Rejected)
        fprintf(stderr,"A total of %d/%d sequences rejected\n",Rejected,*NumEST);
    return HitTable;
}

// ----------------------------
//  Plot Sensor information
// ----------------------------
void SensorEst :: Plot(DNASeq *TheSeq)
{}

// ------------------
//  Post analyse
// ------------------
void SensorEst :: PostAnalyse(Prediction *pred, FILE *MINFO)
{
    int trStart = 0, trEnd = 0, cdsStart = 0, cdsEnd = 0;
    int pprocess = PAR.getI("Est.PostProcess",GetNumber());

    if (!pprocess)
        return;
    if (NumEST == 0)
        return;

    fprintf(MINFO, "#=========================================#\n");
    fprintf(MINFO, "#=             Est evidences             =#\n");
    fprintf(MINFO, "#=========================================#\n");

    qsort((void *)HitTable, NumEST, sizeof(void *), HitsCompareLex);

    // Reset static in EST Support or FEA Support
    if (pprocess == 1)
        ESTSupport(NULL,NULL,100,0,100,0,NULL,0);
    else
        FEASupport(NULL,NULL,100,0,100,0,NULL,0,0);

    for(int i=0; i<pred->nbGene; i++)
    {
        trStart  = pred->vGene[i]->trStart;
        trEnd    = pred->vGene[i]->trEnd;
        cdsStart = pred->vGene[i]->cdsStart;
        cdsEnd   = pred->vGene[i]->cdsEnd;

        if (pprocess == 1)
            // Analyse des ESTs par rapport à la prédiction
            ESTSupport(pred, MINFO, trStart, trEnd, cdsStart, cdsEnd, HitTable, NumEST);
        else
            // Analyse de la prédiction (des features) par rapport aux EST
            FEASupport(pred, MINFO, trStart, trEnd, cdsStart, cdsEnd, HitTable, NumEST, i+1);
    }
}

// -------------------------------------------------------------------------
//  Post analyse (Analyse des ESTs par rapport à la prédiction)
//    Tdebut/fin = debut/fin de transcript
//    debut/fin  = debut/fin de traduit
// -------------------------------------------------------------------------
void SensorEst :: ESTSupport(Prediction *pred, FILE *MINFO, int Tdebut, int Tfin,
                             int debut, int fin,  Hits **HitTable, int Size)
{
    static int EstIndex;

    int supported = 0;
    int CDSsupported = 0;
    unsigned char *Sup;
    Block *ThisBlock;
    int ConsistentEST,i;
    int from = 0, to = 0, ESTEnd = 0;

    if (pred == NULL)
    {
        EstIndex = 0;
        return;
    }

    Sup = new unsigned char[Tfin-Tdebut+1];

    for (i=0; i <= Tfin-Tdebut; i++)
        Sup[i]=0;

    // si la fin du codant n'a pas ete rencontree
    if (fin == -1)
        fin = Tfin;
    if ((debut == -1) || (debut > Tfin))
        debut = Tfin+1;

    //si l'iteration precedente a atteint l'extremite
    if (EstIndex >= Size)
        EstIndex = Max(0,Size-1);

    // on rembobine....
    while ((EstIndex > 0) && (HitTable[EstIndex]->End > Tdebut))
        EstIndex--;

    if (EstIndex >= 0  &&  HitTable[EstIndex]->End < Tdebut)
        EstIndex++;

    while (EstIndex >=0  &&  EstIndex < Size)
    {
        // le dernier transcrit exploitable est passe
        if (HitTable[EstIndex]->Start > Tfin)
            break;

        ConsistentEST = 1;
        ThisBlock = HitTable[EstIndex]->Match;

        while (ThisBlock && ConsistentEST)
        {
            // si on a un gap
            if (ThisBlock->Prev != NULL)
            {
                // intersection
                from = Max(Tdebut,ThisBlock->Prev->End+1);
                to = Min(Tfin,ThisBlock->Start-1);

                for (i = from; i <= to; i++)
                {
                    if (State2Status[pred->GetStateForPos(i+1)] != 1) // 1 = transcribed and spliced
                        ConsistentEST = 0;
                }
            }

            from = Max(Tdebut,ThisBlock->Start);
            to = Min(Tfin,ThisBlock->End);
            ESTEnd = ThisBlock->End;

            for (i = from; i <= to; i++)
            {
                if (State2Status[pred->GetStateForPos(i+1)] < 2) // 2 = transcribed and not spliced
                    ConsistentEST = 0;
            }
            ThisBlock = ThisBlock->Next;
        }
        fprintf(MINFO, "cDNA  %-12s %7d %7d     %4d     %2d introns    ",
                HitTable[EstIndex]->Name,
                HitTable[EstIndex]->Start+1,HitTable[EstIndex]->End+1,
                HitTable[EstIndex]->Length,HitTable[EstIndex]->NGaps);

        if (HitTable[EstIndex]->Rejected)
            fprintf(MINFO, "Filtered ");
        else if (!ConsistentEST)
            fprintf(MINFO, "Inconsistent");
        else
        {
            if ((HitTable[EstIndex]->Start) <= Tdebut && (ESTEnd >= Tfin))
                fprintf(MINFO, "Full Transcript Support");
            else if ((HitTable[EstIndex]->Start) <= debut && (ESTEnd >= fin))
                fprintf(MINFO, "Full Coding Support");
            else
                fprintf(MINFO, "Support");

            ThisBlock = HitTable[EstIndex]->Match;
            while (ThisBlock)
            {
                if (ThisBlock->Prev != NULL)
                {
                    // intersection
                    from = Max(Tdebut,ThisBlock->Prev->End+1);
                    to = Min(Tfin,ThisBlock->Start-1);

                    for (i = from; i <= to; i++)
                        if (!Sup[i-Tdebut])
                        {
                            Sup[i-Tdebut] = 1;
                            supported++;   // no need to test intron/gap consistency (consistent EST)
                            if ((i >=debut) && (i <=fin))
                                CDSsupported++; // !! introns, should not be counted
                        }
                }
                from = Max(Tdebut,ThisBlock->Start);
                to = Min(Tfin,ThisBlock->End);

                for (i= from; i <= to; i++)
                    if (!Sup[i-Tdebut])
                    {
                        Sup[i-Tdebut] = 1;
                        supported++;    // no need to test exon/hit compatibility (consistent EST)
                        if ((i >=debut) && (i <=fin))
                            CDSsupported++;
                    }
                ThisBlock = ThisBlock->Next;
            }
        }
        fprintf(MINFO, "\n");
        EstIndex++;
    }
    if (fin >= debut)
        fprintf(MINFO, "      CDS          %7d %7d    %5d     supported on %d bases\n",
                debut+1,fin+1,fin-debut+1,CDSsupported);
    fprintf(MINFO, "      Gene         %7d %7d    %5d     supported on %d bases\n",
            Tdebut+1,Tfin+1,Tfin-Tdebut+1,supported);
    delete [] Sup;
    return;
}

// -------------------------------------------------------------------------
//  Post analyse (Analyse de la prédiction (features) par rapport aux EST)
//    Tdebut/fin = debut/fin de transcript
//    debut/fin  = debut/fin de traduit
// -------------------------------------------------------------------------
void SensorEst :: FEASupport(Prediction *pred, FILE *MINFO,int Tdebut,int Tfin,
                             int debut,int fin,Hits **HitTable,int Size,int NumGene)
{
    static int EstIndex;

    Block *ThisBlock;
    int ConsistentEST, i, j;
    int from = 0, to = 0, ESTEnd = 0;
    int state;
    std::vector <int> vSupEstI;  // index des transcrits supportant la pred

    if (pred == NULL)
    {
        EstIndex = 0;
        return;
    }


    /***********************************************************************/
    /** Objectif : obtenir un vecteur contenant les index des transcripts **/
    /**            qui supportent la prédiction                           **/
    /***********************************************************************/
    // si la fin du codant n'a pas ete rencontree
    if (fin == -1)
        fin = Tfin;
    if ((debut == -1) || (debut > Tfin))
        debut = Tfin+1;

    // si l'iteration precedente a atteint l'extremite
    if (EstIndex >= Size)
        EstIndex = Max(0, Size-1);

    // on rembobine....
    while ((EstIndex > 0) && (HitTable[EstIndex]->End > Tdebut))
        EstIndex--;

    if (EstIndex >= 0  &&  HitTable[EstIndex]->End < Tdebut)
        EstIndex++;

    while (EstIndex >= 0  &&  EstIndex < Size)
    {
        // le dernier transcrit exploitable est passe
        if (HitTable[EstIndex]->Start > Tfin)
            break;

        ConsistentEST = 1;
        ThisBlock = HitTable[EstIndex]->Match;

        while (ThisBlock && ConsistentEST)
        {
            // si on a un gap
            if (ThisBlock->Prev != NULL)
            {
                // intersection
                from = Max(Tdebut, ThisBlock->Prev->End+1);
                to   = Min(Tfin,   ThisBlock->Start-1);

                for (i = from; i <= to; i++)
                {
                    if (State2Status[pred->GetStateForPos(i+1)] != 1) // 1 = transcrit, épissé
                        ConsistentEST = 0;
                }
            }

            from   = Max(Tdebut, ThisBlock->Start);
            to     = Min(Tfin,   ThisBlock->End);
            ESTEnd = ThisBlock->End;

            for (i = from; i <= to; i++)
            {
                if (State2Status[pred->GetStateForPos(i+1)] < 2) // >= 2 = transcrit non épissé
                    ConsistentEST = 0;
            }
            ThisBlock = ThisBlock->Next;
        }

        // Si EST "coherente"
        if (ConsistentEST)
            vSupEstI.push_back( EstIndex );
        EstIndex++;
    }

    /***********************************************************************/
    /* now, vSupEstI is a set of EST consistent with the analyzed gene     */
    /* Objectif : Analyser chaque feature predite -> supportée ?           */
    /***********************************************************************/
    int  start, end, len;
    char fea[5];
    char strand = pred->vGene[NumGene-1]->vFea[0]->strand;
    Hits **TMPHitTable = new Hits *[NumEST+1]; // need for sorting by % support

    for (i = 0; i < pred->vGene[NumGene-1]->nbFea(); i++)
    {
        state = pred->vGene[NumGene-1]->vFea[i]->state;
        start = pred->vGene[NumGene-1]->vFea[i]->start;
        end   = pred->vGene[NumGene-1]->vFea[i]->end;

        if (end >= Tdebut-1  &&  end <= Tfin+1)
        { // !! assert here ?
            len      = 0;
            int numF = -1;

            if (State2Status[state] >= 2) // >= 2: transcribed unspliced
                numF = pred->vGene[NumGene-1]->vFea[i]->number;

            if (state <= TermR3)
                strcpy(fea, "Exon");
            else if (state == UTR5F  ||  state == UTR5R)
                strcpy(fea, "UTR5");
            else if (state == UTR3F  ||  state == UTR3R)
                strcpy(fea, "UTR3");

            if (numF != -1)
            {
                if (!vSupEstI.empty())
                    // Longueur totale supportée par les transcrits
                    len = LenSup(HitTable, vSupEstI, -1, start, end);

                if (len > 0)
                {
                    fprintf(MINFO, "%s.%d.%d\tEuGene_cDNA\t%s\t%d\t%d\t%d\t%c\t.\t",
                            pred->seqName, (((NumGene-1)*stepid)+1), numF,
                            fea, start, end, len, strand);
                    fprintf(MINFO, "%d\t", (int)((float)len/(end-start+1)*100));

                    for (int k=0;  k<NumEST;  k++)
                        HitTable[k]->Support = 0;

                    // ugly duplicated work here. Would be nice also if EST
                    // consistent with JUST the feature could be analyzed
                    for (j=0; j<(int)vSupEstI.size(); j++)
                    {
                        // Longueur supportée par le transcrit
                        len =  LenSup(HitTable, vSupEstI, j, start, end);
                        HitTable[vSupEstI[j]]->Support = (int)((float)len/(end-start+1)*100);
                    }
                    // On copie la hittable pour trier sur le % supporté
                    for (int k=0; k<NumEST; k++)
                        TMPHitTable[k] = HitTable[k];
                    qsort((void*)TMPHitTable, NumEST, sizeof(void*), HitsCompareSup);

                    // On affiche les ppNumber premiers hits supportant
                    for (j=0; j<NumEST && j<ppNumber && TMPHitTable[j]->Support!=0; j++)
                        fprintf(MINFO, "%s(%d,%d) ", TMPHitTable[j]->Name, TMPHitTable[j]->Support, !TMPHitTable[j]->Rejected);
                    fprintf(MINFO, "\n");
                }
            }
        }
    }

    /***********************************************************************/
    /* Objectif : Analyser la CDS et le gene predit -> supportée ?         */
    /***********************************************************************/
    for (int i=0; i<2; i++)
    {
        if (i==0)
        {
            start = debut;
            end = fin;
            strcpy(fea, "CDS");
        }
        else
        {
            start = Tdebut;
            end = Tfin;
            strcpy(fea, "Gene");
        }

        if (end >= start)
        { // assert ?
            len = 0;
            if (!vSupEstI.empty())
                // Longueur totale supportée par les transcrits
                len = LenSup(HitTable, vSupEstI, -1, start, end);

            if (len > 0)
            {
                fprintf(MINFO, "%s.%d  \tEuGene_cDNA\t%s\t%d\t%d\t%d\t%c\t.\t",
                        pred->seqName, (((NumGene-1)*stepid)+1), fea,
                        start+1, end+1, len, strand);
                fprintf(MINFO, "%d\t", (int)((float)len/(end-start+1)*100));

                for (int k=0;  k<NumEST;  k++)
                    HitTable[k]->Support = 0;

                for (j=0; j<(int)vSupEstI.size(); j++)
                {
                    // Longueur supportée par le transcrit
                    len =  LenSup(HitTable, vSupEstI, j, start, end);
                    HitTable[vSupEstI[j]]->Support = (int)((float)len/(end-start+1)*100);
                }
                // On copie la hittable pour trier sur le % supporté
                for (int k=0; k<NumEST; k++)
                    TMPHitTable[k] = HitTable[k];
                qsort((void*)TMPHitTable, NumEST, sizeof(void*), HitsCompareSup);

                // On affiche les ppNumber premiers hits supportant
                for (j=0; j<NumEST && j<ppNumber && TMPHitTable[j]->Support!=0; j++)
                    fprintf(MINFO, "%s(%d,%d) ", TMPHitTable[j]->Name, TMPHitTable[j]->Support, !TMPHitTable[j]->Rejected);
                fprintf(MINFO, "\n");
            }
        }
    }
    vSupEstI.clear();

    if (TMPHitTable != NULL)
        delete [] TMPHitTable;
    TMPHitTable = NULL;
    return;
}

// -------------------------------------------------------------------------
//  Length supported by EST.
//    If index = -1 -> Lenght supported by ALL EST.
//    Else          -> Lenght supported by ONE EST.
// -------------------------------------------------------------------------
int SensorEst :: LenSup(Hits **HitTable, std::vector<int> vSupEstI,
                        int index, int beg, int end)
{
    int supported = 0;
    unsigned char *Sup;
    Block *ThisBlock;
    int from = 0, to = 0;
    int i;
    int j;

    Sup = new unsigned char[end-beg+1];
    for (i=0; i<=end-beg; i++)
        Sup[i]=0;

    i = (index == -1  ?  0 : index);
    for (; i<(int)vSupEstI.size() || i==index; i++)
    {
        ThisBlock = HitTable[vSupEstI[i]]->Match;

        while (ThisBlock)
        {
            if (ThisBlock->Prev != NULL)
            {
                // intersection
                from = Max(beg, ThisBlock->Prev->End+1);
                to   = Min(end, ThisBlock->Start-1);

                for (j=from; j<=to; j++)
                    if (!Sup[j-beg])
                    {
                        Sup[j-beg] = 1;
                        supported++;
                    }
            }
            from = Max(beg, ThisBlock->Start);
            to   = Min(end, ThisBlock->End);

            for (j=from; j<=to; j++)
                if (!Sup[j-beg])
                {
                    Sup[j-beg] = 1;
                    supported++;
                }
            ThisBlock = ThisBlock->Next;
        }
        if(index != -1)
            break;
    }
    delete [] Sup;
    return supported;
}


void SensorEst :: Print (char name[FILENAME_MAX+1])
{
  FILE *fp;
  strcat (name, ".out");
  if (!(fp = fopen(name, "w"))) {
    fprintf(stderr, "cannot write in %s\n",  name);
    exit(2);
  }

  fprintf(fp, "vPos %d\n",  vPos.size());
  for (int i=0; i< vPos.size();i++ )
  {
    fprintf(fp, "vPos %d\t%d\n",i,  vPos[i]);
  }
  
  fprintf(fp, "vESTMatch %d\n",  vESTMatch.size());
  for (int i=0; i< vESTMatch.size();i++ )
  {
    fprintf(fp, "vESTMatch %d\t\n", vESTMatch[i]);
  }

  fclose(fp);
}
