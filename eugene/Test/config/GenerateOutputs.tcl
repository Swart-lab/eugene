#!/usr/bin/tclsh
#============================================================================
#=             Copyright (c) 2002 by INRA. All rights reserved.             
#=                 Redistribution is not permitted without                  
#=                 the express written permission of INRA.                 
#=                     Mail : tschiex@toulouse.inra.fr                     
#=-------------------------------------------------------------------------
#= File         : EuGeneTk/Test/config/GenerateOutputs.tcl
#= Description  : Generation of the reference files for the test suite
#=                Save the reference output (both stdout and stderr)
#=                for each test
#= Authors      : P.Bardou, S.Foissac, M.J.Cros, A.Moisan, T.Schiex         
#===========================================================================


############################# Environment variables ############################
set action "Generate"

source ./TestProc.tcl
source ./TestVar.tcl
################################################################################


################################################################################
proc AskReplace { file_name } {
	puts "There is a difference between the old reference file $file_name and the new reference file generated.
Replace the old one ? (enter 'Y' for yes else another key for no)" 
}
################################################################################


# Erase all old reference files
puts "Remove all reference files ? (enter 'Y' if you are OK else another key)"
set key [gets stdin]
if {$key=="Y"} {
    foreach f [exec ls ${OUTPUT_DIR}] {
	if {![string match CVS $f]} {
	    exec rm ${OUTPUT_DIR}/$f
	}
    }
    set erase 1
} else {
    set erase 0
}


########################## Init the parameter file ########################
# Convention: if a test requires different parameters values,
#             the initial values are restored after the test
#             BEWARE does not concern sensor use parameters
###########################################################################
# Copy locally the default parameter file
exec cp  ${EUGENE_DIR}/${EUGENE}.par .
# Init parameters values
InitParameterFile ${EUGENE}.par $AllSensorsList


########################################################################
##################        Units tests       ############################
########################################################################
foreach sensor $AllSensorsList {
    # Get stderr and stdout
    if {$sensor != "Est" && $sensor != "Tester"} {
	eval exec $EUGENE_DIR/$EUGENE $OPTIONS(Sensor) \
	    -D Sensor.${sensor}.use=TRUE \
	    $SEQ_DIR/$SEQ(Sensor) 2> tmp%stderr > tmp%stdout
    } else {
	if {$sensor == "Est"} {
	    eval exec $EUGENE_DIR/$EUGENE $OPTIONS(Sensor) \
		-D Sensor.${sensor}.use=TRUE -D Sensor.NG2.use=TRUE \
		$SEQ_DIR/$SEQ(Sensor) 2> tmp%stderr > tmp%stdout
	} else {
	    catch {exec rm test.EuStop.gff}
	    eval exec $EUGENE_DIR/$EUGENE $OPTIONS(Sensor) \
		-D Sensor.${sensor}.use=TRUE $SEQ_DIR/exSeqHom.fasta  \
		2> tmp%stderr > tmp%stdout
	}
    }
    
    # Open files
    set out [open tmp%stdout {RDONLY}]
    set err [open tmp%stderr {RDONLY}]
    set std [open tmp%GenerateOutputs w+]

    # Copy stderr and stdout in the reference file
    set f [read -nonewline $err]
    puts $std $f
    set flux [read -nonewline $out]
    puts $std $flux

    # Close files
    close $out
    close $err
    close $std

    # Remove the 2 first lines related to version number
    RemoveFirstLines tmp%GenerateOutputs
    
    if {$erase == 1 || ![file exists $OUTPUT_DIR/Output_${sensor}]} {
	exec cp tmp%GenerateOutputs $OUTPUT_DIR/Output_${sensor}
    } elseif {[catch {exec diff $OUTPUT_DIR/Output_${sensor} tmp%GenerateOutputs}]} {
	AskReplace $OUTPUT_DIR/Output_${sensor} 
	if {[gets stdin] == "Y"} {
	    exec cp tmp%GenerateOutputs $OUTPUT_DIR/Output_${sensor}
	} else {
	    exec cp tmp%GenerateOutputs $OUTPUT_DIR/Output_${sensor}.new
	}
    }
    
    if {$sensor == "Tester"} {
	if {$erase == 1 || ![file exists $OUTPUT_DIR/Output_${sensor}.gff]} {
	    exec cp test.EuStop.gff $OUTPUT_DIR/Output_${sensor}.gff
	} elseif {[catch {exec diff test.EuStop.gff $OUTPUT_DIR/Output_${sensor}.gff}]} {
	    AskReplace $OUTPUT_DIR/Output_${sensor}.gff
	    if {[gets stdin] == "Y"} {
		eval exec cp test.EuStop.gff $OUTPUT_DIR/Output_${sensor}.gff
	    } else {
		eval exec cp test.EuStop.gff $OUTPUT_DIR/Output_${sensor}.gff.new
	    }
	}
    }

    # Remove all temporary files
    exec rm tmp%GenerateOutputs tmp%stderr tmp%stdout

    puts "Reference files for $sensor unit test created or checked."
}




########################################################################
########################      FUNCTIONAL TESTS      ####################
########################################################################
foreach TEST $FunctionalTestList {
    # Preparation of the parameter file with the correct sensors
    foreach sensor $SensorsList($TEST) \
	{set NewValue${TEST}(Sensor.${sensor}.use) TRUE}
    ModifyParaValue ${EUGENE}.par  NewValue${TEST}

    # Get the sequence length to have only one png file
    set l [GetSeqLength $SEQ_DIR/$SEQ($TEST)]

    # Save output of software and treat them
    eval exec $EUGENE_DIR/$EUGENE $OPTIONS($TEST) -l $l $SEQ_DIR/$SEQ($TEST) \
	2> tmp%stderr > tmp%stdout

    # 1/ image file
    if {$erase == 1 || ![file exists $OUTPUT_DIR/Output_${TEST}.png]} {
	exec cp $IMG($TEST) $OUTPUT_DIR/Output_${TEST}.png
    } elseif {[catch {exec diff $IMG($TEST) $OUTPUT_DIR/Output_${TEST}.png}]} {
	AskReplace $OUTPUT_DIR/Output_${TEST}.png
	if {[gets stdin] == "Y"} {
	    eval exec cp $IMG($TEST) $OUTPUT_DIR/Output_${TEST}.png
	} else {
	    eval exec cp $IMG($TEST) $OUTPUT_DIR/Output_${TEST}.png.new
	}
    }
    # Remove temporary file
    exec rm $IMG($TEST)

    # 2/ reference file in test format
    PrepareReference tmp%stdout tmp%stderr tmp%FunctionalTest

    if {$erase == 1 || ![file exists $OUTPUT_DIR/Output_${TEST}_test]} {
	exec cp tmp%FunctionalTest $OUTPUT_DIR/Output_${TEST}_test
    } elseif {[catch {exec diff tmp%FunctionalTest $OUTPUT_DIR/Output_${TEST}_test}]} {
	AskReplace OUTPUT_DIR/Output_${TEST}_test
	if {[gets stdin] == "Y"} {
	    exec cp tmp%FunctionalTest $OUTPUT_DIR/Output_${TEST}_test
	} else {
	    exec cp tmp%FunctionalTest $OUTPUT_DIR/Output_${TEST}_test.new
	}
    }
    # Remove temporary file
    exec rm tmp%FunctionalTest

    # 3/ reference file in spawn format (no buffered pipes)
    #Open files
    set out [open tmp%stdout {RDONLY}]
    set err [open tmp%stderr {RDONLY}]
    set std [open "tmp%FunctionalTest" w+]
    # Copy stderr and stdout in the reference file
    set f [read -nonewline $err]
    puts $std $f
    set f [read -nonewline $out]
    puts $std $f
    # Close files
    close $out
    close $err
    close $std

    if {$erase == 1 || ![file exists $OUTPUT_DIR/Output_${TEST}]} {
	exec cp tmp%FunctionalTest $OUTPUT_DIR/Output_${TEST}
    } elseif {[catch {exec diff $OUTPUT_DIR/Output_${TEST} tmp%FunctionalTest}]} {
	AskReplace $OUTPUT_DIR/Output_${TEST}
	if {[gets stdin] == "Y"} {
	    exec cp tmp%FunctionalTest $OUTPUT_DIR/Output_${TEST}
	} else {
	    exec cp tmp%FunctionalTest $OUTPUT_DIR/Output_${TEST}.new
	}
    }
    # Remove temporary files
    exec rm tmp%stderr tmp%stdout tmp%FunctionalTest

    # Restore initial parameters values
    InitParameterFile ${EUGENE}.par $AllSensorsList

    puts "Reference files for $TEST functional test created or checked."
}



########################################################################
######################## Sequences base Tests ##########################
########################################################################

# Preparation of the parameter file with the correct sensors
foreach sensor $SensorsList(Araset) \
    {set NewValueAraset(Sensor.${sensor}.use) TRUE}
ModifyParaValue ${EUGENE}.par  NewValueAraset

catch {eval exec $EUGENE_DIR/$EUGENE $SEQ(Araset) > tmp%stdout}

if {$erase == 1 || ![file exists $OUTPUT_DIR/$FILE_REF(Araset)]} {
    exec cp tmp%stdout $OUTPUT_DIR/$FILE_REF(Araset)
} elseif {[catch {exec diff $OUTPUT_DIR/$FILE_REF(Araset) tmp%stdout}]} {
    AskReplace $OUTPUT_DIR/$FILE_REF(Araset)
    if {[gets stdin] == "Y"} {
	exec cp tmp%stdout $OUTPUT_DIR/$FILE_REF(Araset)
    } else {
	exec cp tmp%stdout $OUTPUT_DIR/${FILE_REF(Araset)}.new
    }
}

# remove temporary file	
exec rm tmp%stdout

# Restore initial parameters values
InitParameterFile ${EUGENE}.par $AllSensorsList

puts "Reference files for Araset test created or checked."



########################################################################
##################### Parameters optimization ##########################
########################################################################

catch {eval exec $EUGENE_DIR/$EUGENE test -D ParaOptimization.Use=TRUE > tmp%stdout}

if {$erase == 1 || ![file exists $OUTPUT_DIR/$FILE_REF(Optimization)]} {
    exec cp tmp%stdout $OUTPUT_DIR/$FILE_REF(Optimization)
} elseif {[catch {exec diff $OUTPUT_DIR/$FILE_REF(Optimization) tmp%stdout}]} {
    AskReplace $OUTPUT_DIR/$FILE_REF(Optimization)
    if {[gets stdin] == "Y"} {
	exec cp tmp%stdout $OUTPUT_DIR/$FILE_REF(Optimization)
    } else {
	exec cp tmp%stdout $OUTPUT_DIR/$FILE_REF(Optimization).new
    }
}
# remove temporary file	
exec rm tmp%stdout

puts "Reference files for Optimization test created or checked."


# Indicate the end of the reference files generation
puts "Reference files generated in the $OUTPUT_DIR directory."

catch {eval exec rm ./$EUGENE.par}

