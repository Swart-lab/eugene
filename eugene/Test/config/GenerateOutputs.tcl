#!/usr/bin/tcl
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
	exec rm ${OUTPUT_DIR}/$f
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
    if {$sensor != "Est"} {
	eval exec $EUGENE_DIR/$EUGENE $SEQ_DIR/$SEQ(Sensor) $OPTIONS(Sensor) \
	    -D Sensor.${sensor}.use=TRUE 2> stderr > stdout
    } else {
	eval exec $EUGENE_DIR/$EUGENE $SEQ_DIR/$SEQ(Sensor) $OPTIONS(Sensor) \
	    -D Sensor.${sensor}.use=TRUE -D Sensor.NG2.use=TRUE 2> stderr > stdout
    }

# Open files
    set out [open stdout {RDONLY}]
    set err [open stderr {RDONLY}]
    set std [open GenerateOutputs.tmp w+]

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
    RemoveFirstLines GenerateOutputs.tmp
    
    if {$erase == 1 || ![file exists $OUTPUT_DIR/Output_${sensor}]} {
	exec cp GenerateOutputs.tmp $OUTPUT_DIR/Output_${sensor}
    } elseif {[catch {exec diff $OUTPUT_DIR/Output_${sensor} GenerateOutputs.tmp}]} {
	AskReplace $OUTPUT_DIR/Output_${sensor} 
	if {[gets stdin] == "Y"} {
	    exec cp GenerateOutputs.tmp $OUTPUT_DIR/Output_${sensor}
	} else {
	    exec cp GenerateOutputs.tmp $OUTPUT_DIR/Output_${sensor}.new
	}
    }
	
# Remove all temporary files
    exec rm GenerateOutputs.tmp stderr stdout
}




########################################################################
########################      FUNCTIONAL TESTS      ####################
########################################################################
foreach TEST $FunctionalTestList {
    # Preparation of the parameter file with the correct sensors
    foreach sensor $SensorsList($TEST) {set NewValue${TEST}(Sensor.${sensor}.use) TRUE}
    ModifyParaValue ${EUGENE}.par  NewValue${TEST}

    # Get the sequence length
    set l [GetSeqLength $SEQ_DIR/$SEQ($TEST)]

    # Save output of software and treat them
    eval exec $EUGENE_DIR/$EUGENE $SEQ_DIR/$SEQ($TEST) $OPTIONS($TEST) -l $l 2> stderr > stdout

    # 1/ image file
    if {$erase == 1 || ![file exists $OUTPUT_DIR/Output_${TEST}.000.png]} {
	exec cp $IMG($TEST) $OUTPUT_DIR/Output_${TEST}.000.png
    } elseif {[catch {exec diff $IMG($TEST) $OUTPUT_DIR/Output_${TEST}.000.png}]} {
	AskReplace $OUTPUT_DIR/Output_${TEST}.000.png
	if {[gets stdin] == "Y"} {
	    eval exec cp $IMG($TEST) $OUTPUT_DIR/Output_${TEST}.000.png
	} else {
	    eval exec cp $IMG($TEST) $OUTPUT_DIR/Output_${TEST}.000.png.new
	}
    }
    # Remove temporary file
    exec rm $IMG($TEST)

    # 2/ reference file in test format
    PrepareReference stdout stderr FunctionalTest.tmp 

    if {$erase == 1 || ![file exists $OUTPUT_DIR/Output_${TEST}]} {
	exec cp FunctionalTest.tmp $OUTPUT_DIR/Output_${TEST}
    } elseif {[catch {exec diff FunctionalTest.tmp $OUTPUT_DIR/Output_${TEST}}]} {
	AskReplace OUTPUT_DIR/Output_${TEST}
	if {[gets stdin] == "Y"} {
	    exec cp FunctionalTest.tmp $OUTPUT_DIR/Output_${TEST}
	} else {
	    exec cp FunctionalTest.tmp $OUTPUT_DIR/Output_${TEST}.new
	}
    }
    # Remove temporary file
    exec rm FunctionalTest.tmp

    # 3/ reference file in spawn format (no buffered pipes)
    #Open files
    set out [open stdout {RDONLY}]
    set err [open stderr {RDONLY}]
    set std [open "FunctionalTest.tmp" w+]
    # Copy stderr and stdout in the reference file
    set f [read -nonewline $err]
    puts $std $f
    set f [read -nonewline $out]
    puts $std $f
    # Close files
    close $out
    close $err
    close $std

    if {$erase == 1 || ![file exists $OUTPUT_DIR/Output_${TEST}_std]} {
	exec cp FunctionalTest.tmp $OUTPUT_DIR/Output_${TEST}_std
    } elseif {[catch {exec diff $OUTPUT_DIR/Output_${TEST}_std FunctionalTest.tmp}]} {
	AskReplace $OUTPUT_DIR/Output_${TEST}_std
	if {[gets stdin] == "Y"} {
	    exec cp FunctionalTest.tmp $OUTPUT_DIR/Output_${TEST}_std
	} else {
	    exec cp FunctionalTest.tmp $OUTPUT_DIR/Output_${TEST}_std.new
	}
    }
    # Remove temporary files
    exec rm stderr stdout FunctionalTest.tmp

    # Restore initial parameters values
    InitParameterFile ${EUGENE}.par $AllSensorsList
}



########################################################################
######################## Sequences base Tests ##########################
########################################################################

# Preparation of the parameter file with the correct sensors
foreach sensor $SensorsList(Araset) {set NewValueAraset(Sensor.${sensor}.use) TRUE}
ModifyParaValue ${EUGENE}.par  NewValueAraset

catch {eval exec $EUGENE_DIR/$EUGENE $SEQ(Araset) > stdout}

if {$erase == 1 || ![file exists $OUTPUT_DIR/Output_Araset]} {
    exec cp stdout $OUTPUT_DIR/Output_Araset
} elseif {[catch {exec diff $OUTPUT_DIR/Output_Araset stdout}]} {
    AskReplace $OUTPUT_DIR/Output_Araset
    if {[gets stdin] == "Y"} {
	exec cp stdout $OUTPUT_DIR/Output_Araset
    } else {
	exec cp stdout $OUTPUT_DIR/Output_Araset.new
    }
}

# remove temporary file	
exec rm stdout

# Restore initial parameters values
InitParameterFile ${EUGENE}.par $AllSensorsList




# Indicate the end of the reference files generation
puts "Reference files generated in the $OUTPUT_DIR directory."

catch {eval exec rm ./$EUGENE.par}