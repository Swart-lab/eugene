#============================================================================
#=             Copyright (c) 2002 by INRA. All rights reserved.             
#=                 Redistribution is not permitted without                  
#=                 the express written permission of INRA.                 
#=                     Mail : tschie@toulouse.inra.fr                     
#=-------------------------------------------------------------------------
#= File         : EuGeneTk/Test/config/TestVar.exp
#= Description  : Environment variables definition
#= Authors      : P.Bardou, S.Foissac, M.J.Cros, A.Moisan, T.Schiex         
#===========================================================================


############################# Environment variables ###########################
set EUGENE EuGeneTest

if {$action=="Test"} {
    set EUGENE_DIR .
    set OUTPUT_DIR ../Test/Outputs
    set SEQ_DIR ../Test/Sequences
    set TRACE_DIR ../Test/TestTrace
} elseif {$action=="Generate"} {
    set EUGENE_DIR ../../EuGene
    set OUTPUT_DIR ../Outputs
    set SEQ_DIR ../Sequences
    set TRACE_DIR ../TestTrace
}

############################################################################
set AllSensorsList {ATGpr BlastX Est EuStop FrameShift GFF GSplicer Homology IfElse \
		     MarkovConst MarkovIMM MarkovProt NG2 NStart PatConst Repeat  \
		     Riken SPred SpliceWAM \
		     StartWAM Transcript User GCPlot Plotter Tester}


############################################################################
##################### Units tests variables ################################
set SEQ(Sensor) {seq14ac002535g4g5.tfa}
set OPTIONS(Sensor) "-pd"


############################################################################# 
set FunctionalTestList {SeqAra SeqDoc SeqHom SeqRest SeqBig}

##################### SeqAra test variables #################################
set SensorsList(SeqAra) {MarkovIMM MarkovConst EuStop NStart IfElse GSplicer Est BlastX}
set SEQ(SeqAra) {seq25ab005234g10g11.tfa}
set IMG(SeqAra) {seq25ab005234g10g11.test.000.png}
set FILE_REF(SeqAra) Output_SeqAra
set OPTIONS(SeqAra) "-gtest -E -B"

##################### SeqDoc test variables ############################
set SensorsList(SeqDoc) {MarkovIMM MarkovConst EuStop NStart IfElse Transcript \
			     FrameShift User}
set SEQ(SeqDoc) {SYNO_ARATH.fasta}
set IMG(SeqDoc) {SYNO_ARATH.test.000.png}
set FILE_REF(SeqDoc) Output_SeqDoc
set OPTIONS(SeqDoc) "-gtest"

##################### SeqHom test variables ################################
set SensorsList(SeqHom) {MarkovProt MarkovConst EuStop StartWAM SpliceWAM Homology \
			     Transcript User GFF}
set SEQ(SeqHom) {exSeqHom.fasta}
set IMG(SeqHom) {exSeqHom.test.000.png}
set FILE_REF(SeqHom) Output_SeqHom
set OPTIONS(SeqHom) "-gtest"
set NewValueSeqHom(Transcript.Start*) 1e-9
set NewValueSeqHom(Transcript.Stop*)  1e-9
set NewValueSeqHom(EuStop.stopP*)     3.1439
set NewValueSeqHom(BlastX.level1*)    0.0

##################### SeqRest test variables #################################
set SensorsList(SeqRest) {MarkovIMM EuStop PatConst IfElse Riken Repeat}
set SEQ(SeqRest) {seq25ab005234g10g11.tfa}
set IMG(SeqRest) {seq25ab005234g10g11.test.000.png}
set FILE_REF(SeqRest) Output_SeqRest
set OPTIONS(SeqRest) "-gtest"

##################### SeqBig test variables ################################
set SensorsList(SeqBig) {MarkovIMM MarkovConst EuStop NStart IfElse GSplicer Est BlastX}
set SEQ(SeqBig) {chr1_2002.tfa_22540000-23040000.tfa}
set IMG(SeqBig) {chr1_2002.tfa_22540000-23040000.test.000.png}
set FILE_REF(SeqBig) Output_SeqBig
set OPTIONS(SeqBig) "-gtest"


############################################################################
##################### Araset test variables ################################
set SEQBASE_DIR "/Annotation/Arabidopsis/araset/Genes"
set SensorsList(Araset) {MarkovIMM MarkovConst EuStop NStart IfElse GSplicer Est BlastX}
set SEQ(Araset) {}
set seq [eval exec ls $SEQBASE_DIR]
foreach f $seq {
    if { [string match *.tfa $f] } {
	set ff ${SEQBASE_DIR}/$f
	lappend SEQ(Araset) $ff
    }
}
set FILE_REF(Araset) Output_Araset



############################################################################
##################### Parameters optimization ##############################
set FILE_REF(Optimization) Output_Optimization



