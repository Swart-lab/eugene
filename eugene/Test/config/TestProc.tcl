#============================================================================
#=             Copyright (c) 2002 by INRA. All rights reserved.             
#=                 Redistribution is not permitted without                  
#=                 the express written permission of INRA.                 
#=                     Mail : tschiex@toulouse.inra.fr                     
#=-------------------------------------------------------------------------
#= File         : EuGeneTk/Test/Reference/TestProc.tcl
#= Description  : Procedures for test of EuGene sofware
#= Authors      : P.Bardou, S.Foissac, M.J.Cros, A.Moisan, T.Schiex         
#===========================================================================



############################################################################
# Procedure ModifyParaValue
# Description : in the parameter file specified by the 1st argument, 
#               modify the value of the parameters given as index 
#               of the 2sd argument (array)
#               which contains the new values of parameters
# Example     : tcl> set V(EuGene.minEx) 4
#               tcl> set V(Output.graph) TRUE
#               tcl> ModifyParaValue EuGeneAS.par V
#               After execution, in the file EuGeneAS.par 
#               the line 'EuGene.minEx   3' will be ' EuGene.minEx   4' 
#               the line 'Output.graph   FALSE' will be 'Output.graph  TRUE'   
# Note that if there is no ambiguity, it is sufficient to give 
# the end of parameter name
# In the example, set V(minEx) 4 would produces the same result 
# if there is not an other 'minEx' in the file.
# BEWARE in case of parameter name included in an other: to distinguish
# the shortest it is necessary to add a space at the end of its name
############################################################################
proc ModifyParaValue {FileName NewValues} {
# To pass an array as 2nd argument
    upvar $NewValues V

    set f [open $FileName r]
    set new_content ""

    foreach line [split [read $f] \n] {
	if {[string length $line] != 0} {
	    set newline ""
	    foreach para [array names V] {
		set Pos [string first ${para} $line]
		if {$Pos != -1} {
		    set newline [string range $line \
				     0 [expr $Pos - 1 + [string length $para]]]
		    set newline "$newline $V($para)"
		}
	    }
	    if {[string length $newline] == 0} {
		set newline $line
	    }
	    set new_content "$new_content$newline\n"
	}
    }
    close $f

# Write the new parameter file
    set f [open $FileName w]
    puts -nonewline $f $new_content
    close $f
}


##############################################################################
# Procedure   : InitParameterFile
# Description : Update values of parameters in a parameter file
#               Arguments: FileName = name of the parameter file
#                          SensorsList = list of sensors that will be specified 
#                                        not to use
# BEWARE      : Priority of sensors are not updated
###############################################################################
proc InitParameterFile {FileName SensorsList} {

# when name of parameter is include in an other
# (example SpliceConst.accP is included in SpliceConst.accPNo)
# it is necessary to add a space at the name of the sensor
# (example SpliceConst.accP$space)
set space " "

# Set the parameter to wanted values
set NewValue1(EuGene.minEx)	3
set NewValue1(EuGene.minIn)	50
set NewValue1(EuGene.minSg)	150
set NewValue1(EuGene.minFlow)	50
set NewValue1(EuGene.minConv)	0
set NewValue1(EuGene.minDiv)	100
set NewValue1(EuGene.min5Prime)	20
set NewValue1(EuGene.min3Prime)	20
##### Priors #####
set NewValue1(EuGene.ExonPrior)	0.33	
set NewValue1(EuGene.IntronPrior)	0.17	
set NewValue1(EuGene.InterPrior)	0.4
set NewValue1(EuGene.FivePrimePrior)	0.03
set NewValue1(EuGene.ThreePrimePrior)	0.07
##### States penalties #####
set NewValue1(EuGene.transCodant)	1.0
set NewValue1(EuGene.transIntron)	1.0
set NewValue1(EuGene.transUTR5)	0.999
set NewValue1(EuGene.transUTR3)	0.999
set NewValue1(EuGene.transInter)	1.0
##### Output control ######
set NewValue1(Output.graph)		FALSE	
set NewValue1(Output.resx)		900
set NewValue1(Output.resy)		400
set NewValue1(Output.glen)		-1
set NewValue1(Output.golap)		-1
set NewValue1(Output.gfrom)		-1
set NewValue1(Output.gto)		-1
set NewValue1(Output.window)		48
set NewValue1(Output.format)		l
set NewValue1(Output.offset)		0
set NewValue1(Output.normopt)		1
set NewValue1(Output.Prefix)		./
##### Transcript parameters #####
set NewValue1(Transcript.Start)	4.155
set NewValue1(Transcript.Stop)		4.155
##### EuStop parameters #####
set NewValue1(EuStop.stopP)            4.155
##### FrameShift parameters #####
set NewValue1(FrameShift.Ins)	999.0
set NewValue1(FrameShift.Del)	999.0
##### NetStart parameters #####
set NewValue1(NStart.startP)	0.052
set NewValue1(NStart.startB)	0.308
##### IfElse #####
set NewValue1(IfElse.SensorIf)		NG2
set NewValue1(IfElse.SensorElse)	SPred
# NetGene2 parameters #####
set NewValue1(NG2.accP\[0\])     0.903
set NewValue1(NG2.accB\[0\])     5.585
set NewValue1(NG2.donP\[0\])     0.980
set NewValue1(NG2.donB\[0\])     27.670
set NewValue1(NG2.accP\[1\])	0.903
set NewValue1(NG2.accB\[1\])	5.585
set NewValue1(NG2.donP\[1\])	0.980
set NewValue1(NG2.donB\[1\])	27.670
##### SplicePredictor parameters #####
set NewValue1(SPred.accP\[0\])   0.987
set NewValue1(SPred.accB\[0\])  3.850
set NewValue1(SPred.donP\[0\])   0.929
set NewValue1(SPred.donB\[0\])   10.800
set NewValue1(SPred.accP\[1\])	0.987
set NewValue1(SPred.accB\[1\])	3.850
set NewValue1(SPred.donP\[1\])	0.929
set NewValue1(SPred.donB\[1\])	10.800
##### Interpolated Markov Models parameters #####
set NewValue1(MarkovIMM.matname\[0\])	Ara2UTR.mat
set NewValue1(MarkovIMM.minGC\[0\])	0
set NewValue1(MarkovIMM.maxGC\[0\])	100
##### Markov proteic model parameters #####
set NewValue1(MarkovProt.matname\[0\])	swissprot.maxorder2.bin
set NewValue1(MarkovProt.minGC\[0\])	0
set NewValue1(MarkovProt.maxGC\[0\])	100
set NewValue1(MarkovProt.maxorder)    2
set NewValue1(MarkovProt.order)        2
##### Est sensor parameters #####
set NewValue1(Est.PostProcess)	FALSE
set NewValue1(Est.estP)	-0.4
set NewValue1(Est.estM)	6
set NewValue1(Est.utrP)	0.35
set NewValue1(Est.utrM)	5
##### Proteic similarity sensor parameters #####
set NewValue1(BlastX.PostProcess) FALSE
set NewValue1(BlastX.levels)	0
set NewValue1(BlastX.level0)	0.2
set NewValue1(BlastX.level1)	0.05
set NewValue1(BlastX.level2)	0.0
set NewValue1(BlastX.level3)	0.0
set NewValue1(BlastX.level4)	0.0
set NewValue1(BlastX.level5)	0.0
set NewValue1(BlastX.level6)	0.0
set NewValue1(BlastX.level7)	0.0
set NewValue1(BlastX.level8)	0.0
set NewValue1(BlastX.level9)	0.0
set NewValue1(BlastX.blastxM)	10	
##### Repeat sensor parameters #####
set NewValue1(Repeat.UTRPenalty)	0.0
set NewValue1(Repeat.IntronPenalty)	0.1
set NewValue1(Repeat.ExonPenalty)	1.0
##### Homology Sensor parameters #####
set NewValue1(Homology.TblastxP) 	0
set NewValue1(Homology.TblastxB) 	0.0595
set NewValue1(Homology.protmatname	BLOSUM80
##### Sensors with uniform penalties #####
set NewValue1(SpliceConst.accP$space)	1e-9
set NewValue1(SpliceConst.accPNo)       0
set NewValue1(SpliceConst.donP$space)    1e-9
set NewValue1(SpliceConst.donPNo)       0
set NewValue1(StartConst.startP$space)	2.897949
set NewValue1(StartConst.startPNo)	0
set NewValue1(MarkovConst.val)		0.25
set NewValue1(MarkovConst.minGC\[0\])	0
set NewValue1(MarkovConst.maxGC\[0\])	100
#
##### Sensor SpliceWAM #####
#
set NewValue1(SpliceWAM.MarkovianOrder)	1
set NewValue1(SpliceWAM.donmodelfilename)	WAM/WAM.ARA.DON.L9
set NewValue1(SpliceWAM.NbNtBeforeGT)	3
set NewValue1(SpliceWAM.NbNtAfterGT)	4
set NewValue1(SpliceWAM.DonScaleCoef)	2.9004
set NewValue1(SpliceWAM.DonScalePenalty)	-7.5877
set NewValue1(SpliceWAM.accmodelfilename)	WAM/WAM.ARA.ACC.L7
set NewValue1(SpliceWAM.NbNtBeforeAG)		2
set NewValue1(SpliceWAM.NbNtAfterAG)		1
set NewValue1(SpliceWAM.AccScaleCoef)		2.9004
set NewValue1(SpliceWAM.AccScalePenalty)		-7.5877
#
##### Sensor StartWAM #####
#
set NewValue1(StartWAM.modelfilename)	WAM/WAM.ARA.START9
set NewValue1(StartWAM.NbNtBeforeATG)	3
set NewValue1(StartWAM.NbNtAfterATG)	3
set NewValue1(StartWAM.MarkovianOrder)		1
set NewValue1(StartWAM.ScaleCoef)		0.1594
set NewValue1(StartWAM.ScalePenalty)		-3.1439


# No sensor
foreach sensor $SensorsList {set NewValue1(Sensor.${sensor}.use) FALSE}

ModifyParaValue $FileName  NewValue1
}




#############################################################################
# Procedure   : BackQuoteLine
# Description : Add a '\\' before '+' '(' ')' '|' '[' ']'
#               Argument : OldLine = string to consider
# Evaluation  : modified string 
#############################################################################
proc BackQuoteLine {OldLine} {
    set l $OldLine
    foreach sign { + ( ) | [ ]} {
	set NewLine ""
	set sign_pos [string first $sign $l]
	while { $sign_pos != -1} {
	    set NewLine "$NewLine[string range $l 0 [expr $sign_pos - 1]]\\"
	    set NewLine "$NewLine$sign"
	    set l [string range $l [expr $sign_pos + 1] [string length $l]]
	    set sign_pos [string first $sign $l]
	}
	set NewLine "$NewLine$l"
	set l $NewLine
    }
    return $NewLine
}


#############################################################################
# Procedure   : PrepareReference
# Description : Modify an output file in a reference file for test:
#                 - remove the 2 first lines related to version number
#                 - move the end of stderr bad placed at the end of the file
#                 - back quote '+', '(', ')', '|'
#               Argument : FileName = name of an output file obtained
#                          with a command EuGeneAS -g ... >& FileName
#                          (write first all stderr, then all stdout)
# BEWARE      : MUST be an output with graph (-g command argument)
#############################################################################
proc PrepareReference {stdout stderr FileName} {

# Remove the first lines related to version number
    RemoveFirstLines $stderr

# Copy stderr and stdout in the good order in a FileName

    set out [open $stdout {RDONLY}]
    set err [open $stderr {RDONLY}]
    set f [open $FileName w]
    set content ""

    set newline [gets $err]
    while { ![string match "*Seq         Type    S*" $newline] } {
	set content "${content}${newline}\n"
	set newline [gets $err]
    }
    set content "${content}${newline}\n"
    set newline [gets $err]
    set content "${content}${newline}\n"

    set newline [read -nonewline $out]
    set content "${content}${newline}\n"

    set newline [read -nonewline $err]
    set content "${content}${newline}"

    puts $f $content

    close $f
    close $out
    close $err

# Place a '.' at the end of line
# Back quote '+', '(', ')', '|'
    set f [open $FileName r]
    set new_content ""
    foreach line [split [read -nonewline $f] \n] {
	set line [BackQuoteLine $line]
	set new_content "${new_content}${line}.\n"
    }
    close $f
    set f [open $FileName w]
    puts $f $new_content
    close $f
}


#############################################################################
# Procedure   : RemoveFirstLines
# Description : Remove the 2 first lines in the file given in argument
#############################################################################
proc RemoveFirstLines {file_name} {
   exec cp $file_name RemoveFirstLines.tmp
   exec tail +3 RemoveFirstLines.tmp > $file_name
   exec rm RemoveFirstLines.tmp
}


#############################################################################
# Procedure   : GetSeqLength
# Description : returns the number of caracters of the file given in argument
#############################################################################
proc GetSeqLength {file_name} {
    set l [exec wc -c $file_name]
    set l [string trim $l]
    set l [lindex [split $l] 0]
    return $l
}

