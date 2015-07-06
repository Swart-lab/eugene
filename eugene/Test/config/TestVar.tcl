# ------------------------------------------------------------------
# Copyright (C) 2004 INRA <eugene@ossau.toulouse.inra.fr>
#
# This program is open source; you can redistribute it and/or modify
# it under the terms of the Artistic License (see LICENSE file).
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
#
# You should have received a copy of Artistic License along with
# this program; if not, please see http://www.opensource.org
#
# $Id$
# ------------------------------------------------------------------
# File:     TestVar.tcl
# Contents: Environment variables definition
# ------------------------------------------------------------------


############################# Environment variables ###########################
set EUGENE_TEST_PAR      eugeneTest.par
set EUGENE_TEST_PROK_PAR eugeneTestProk.par
set EUGENE src/eugene

if {$action=="Test"} {
# ask tests in the root directory
    set EUGENE_DIR [exec pwd]
    set OUTPUT_DIR Test/Outputs
    set SEQ_DIR Test/Sequences
    set TRACE_DIR Test/TestTrace
} elseif {$action=="Generate"} {
# ask generation of ouput files in Test/config
    set EUGENE_DIR [exec pwd]/../..
    set OUTPUT_DIR ../Outputs
    set SEQ_DIR ../Sequences
    set TRACE_DIR ../TestTrace
}

############################################################################
set AllSensorsList {AnnotaStruct BlastX Est EuStop FrameShift GCPlot GFF GSplicer \
		    Homology IfElse MarkovConst MarkovIMM MarkovProt NG2 NStart PatConst  \
		    PepSignal Plotter Repeat Riken SMachine SPred  \
		    Tester Transcript NcRNA ProStart}

############################################################################
# At the moment, NcRNA  is not in this list because the gff3 format for this 
# sensor is already tested in unit tests
set AllGff3Sensors {AnnotaStruct BlastX Est GFF GSplicer NG2 NStart \
		    PepSignal Repeat SMachine SPred Homology}
############################################################################
##################### Units tests variables ################################
set SEQ(Sensor) {seq14ac002535g4g5}
set OPTIONS(Sensor) "-pd"

############################################################################# 
set FunctionalTestList {SeqAra SeqDoc SeqAlt SeqRest SeqNcRNA ProOverlapGene EstNonCanSite}

##################### ProOverlapGene test variables #########################
set SensorsList(ProOverlapGene) {MarkovIMM ProStart EuStop Transcript BlastX}
set SEQ(ProOverlapGene) {SMc.1541000-1552500.fasta}
set IMG(ProOverlapGene) {SMc.1541000-1552500.000.png}
set GFF3(ProOverlapGene) {SMc.1541000-1552500.gff3}
set FILE_REF(ProOverlapGene) Output_ProOverlapGene
set OPTIONS(ProOverlapGene) "-pog -g"
set NewValueProOverlapGene(MarkovIMM.matname\[0\]) Sm.mat	
set NewValueProOverlapGene(EuStop.stopP*) 4.6644
set NewValueProOverlapGene(Transcript.Start*)	   4.6644
set NewValueProOverlapGene(Transcript.Stop*)	   4.6644
set NewValueProOverlapGene(Transcript.StartNpc*)   4.6644
set NewValueProOverlapGene(Transcript.StopNpc*)    4.6644
set NewValueProOverlapGene(Transcript.AffectedStrand) 0
set NewValueProOverlapGene(BlastX.levels)	   01
set NewValueProOverlapGene(BlastX.level0*)	0.7756
set NewValueProOverlapGene(BlastX.level1*)	0.0506
set NewValueProOverlapGene(ProStart.alpha*)	0.8733
set NewValueProOverlapGene(ProStart.beta*)	7.7582
set NewValueProOverlapGene(Sensor.MarkovConst.use) 0

##################### SeqAra test variables #################################
set SensorsList(SeqAra) {MarkovIMM MarkovConst EuStop NStart IfElse GSplicer Est BlastX}
set SEQ(SeqAra) {seq25ab005234g10g11.tfa}
set IMG(SeqAra) {seq25ab005234g10g11.000.png}
set FILE_REF(SeqAra) Output_SeqAra
set OPTIONS(SeqAra) "-po -g"

##################### SeqDoc test variables ############################
set SensorsList(SeqDoc) {MarkovIMM MarkovConst EuStop NStart IfElse Transcript \
		     FrameShift }
set SEQ(SeqDoc) {SYNO_ARATH.fasta}
set IMG(SeqDoc) {SYNO_ARATH.000.png}
set FILE_REF(SeqDoc) Output_SeqDoc
set OPTIONS(SeqDoc) "-po -g"

##################### SeqHom test variables ################################
set SensorsList(SeqHom) {MarkovProt MarkovConst EuStop StartWAM SpliceWAM Homology \
			     Transcript GFF}
set SEQ(SeqHom) {exSeqHom.fasta}
set IMG(SeqHom) {exSeqHom.000.png}
set FILE_REF(SeqHom) Output_SeqHom
set OPTIONS(SeqHom) "-po -g"
set NewValueSeqHom(Transcript.Start*) 1e-9
set NewValueSeqHom(Transcript.Stop*)  1e-9
set NewValueSeqHom(Transcript.AffectedStrand) 0
set NewValueSeqHom(EuStop.stopP*)     3.1439

##################### SeqRest test variables #################################
set SensorsList(SeqRest) {MarkovIMM EuStop PatConst IfElse Riken Repeat}
set SEQ(SeqRest) {seq25ab005234g10g11.tfa}
set IMG(SeqRest) {seq25ab005234g10g11.000.png}
set FILE_REF(SeqRest) Output_SeqRest
set OPTIONS(SeqRest) "-po -g"

##################### EstNonCanSite test variables #########################
set SensorsList(EstNonCanSite) {EuStop NStart SPred Transcript Est}
set SEQ(EstNonCanSite) {seq_noncansplicesites.tfa}
set IMG(EstNonCanSite) {seq_noncansplicesites.000.png}
set FILE_REF(EstNonCanSite) Output_EstNonCanSite
set OPTIONS(EstNonCanSite) "-po -g -DEuGene.NonCanDon=TT -DEuGene.NonCanAcc=GG "
set NewValueEstNonCanSite(Est.SpliceNonCanP\[0\]) 1.0


############################################################################
set ArabidopsisTestList {Araset Big}

##################### Araset test variables ################################
set SEQBASE_DIR "/Annotation/Arabidopsis/araset/Genes"
#"/home/bardou/araset/Genes"
set SensorsList(Araset) {MarkovIMM MarkovConst EuStop NStart IfElse GSplicer Est BlastX Transcript}
set SEQ(Araset) {}
set seq [eval exec ls $SEQBASE_DIR]
foreach f $seq {
    if { [string match *.tfa $f] } {
	set ff ${SEQBASE_DIR}/$f
	lappend SEQ(Araset) $ff
    }
}
set FILE_REF(Araset) Output_Araset
set OPTIONS(Araset) "-po"
set FILE_REF(ArasetSpSn) Output_Araset_SpSn 
set FILE_COORD(ArasetSpSn) "/Annotation/Arabidopsis/araset/Truth/araset.coord"
set PRG_EVAL_PRED $EUGENE_DIR/Procedures/Eval/egn_evalpred.pl

##################### Big test variables ################################
set SensorsList(Big) {MarkovIMM  Transcript MarkovConst EuStop NStart  GSplicer Est BlastX}
set SEQ(Big) $SEQ_DIR/chr1_2002.tfa_22540000-23040000.tfa
set FILE_REF(Big) Output_Big
set OPTIONS(Big) "-po"

##################### SeqAlt test variables #################################
set SensorsList(SeqAlt) {MarkovIMM MarkovConst EuStop NStart IfElse Est NG2 SPred}
set SEQ(SeqAlt) {At5g18830.fasta.genomicAJ011613.fasta}
set IMG(SeqAlt) {At5g18830.fasta.genomicAJ011613.000.png}
set FILE_REF(SeqAlt) Output_SeqAlt
set OPTIONS(SeqAlt) "-po -g -a"

##################### SeqNcRNA test variables ###############################
set SensorsList(SeqNcRNA) {MarkovIMM MarkovConst NStart Transcript NcRNA EuStop IfElse}
set SEQ(SeqNcRNA) {seq14ac002535g4g5.tfa}
set IMG(SeqNcRNA) {seq14ac002535g4g5.000.png}
set FILE_REF(SeqNcRNA) Output_SeqNcRNA
set NewValueSeqNcRNA(NcRNA.NpcRna*)     0.1
set NewValueSeqNcRNA(NcRNA.TStartNpc*) 100
set NewValueSeqNcRNA(NcRNA.TStopNpc*) 100
set OPTIONS(SeqNcRNA) "-po -g -a"
######################################################################################### Parameters optimization ##############################
set FILE_REF(Optimization) Output_Optimization


