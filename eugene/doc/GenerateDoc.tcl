#!/usr/bin/tcl
#============================================================================
#=             Copyright (c) 2002 by INRA. All rights reserved.             
#=                 Redistribution is not permitted without                  
#=                 the express written permission of INRA.                 
#=                     Mail : tschiex@toulouse.inra.fr                     
#=-------------------------------------------------------------------------
#= File         : EuGeneTk/Documentation/GenerateDoc.tcl               
#= Description  : Generation of the EuGene documentation
#= Authors      : P.Bardou, S.Foissac, M.J.Cros, A.Moisan, T.Schiex         
#===========================================================================



#===========================================================================
# definition of variables related to I/O

set SEQ Sequences/SYNO_ARATH.fasta
set EUGENEAS_DIR ../EuGene/
set EUGENEAS EuGeneAS 
exec cp ../EuGene/EuGeneAS.par .

set FIC_TEX "Doc.tex"
set FIC_TEX_TMP "EuGeneDoc"
set FIC_TMP "tmp"
set CMDFLAGS_INDEX "CmdFlags"

set Cmd_end " >& $FIC_TMP"
set Flag(1) EXECUTION_TRACE1; set Cmd_begin(1) ""; set Cmd(1) "$EUGENEAS $SEQ -s"
set Flag(2) EXECUTION_TRACE2; set Cmd_begin(2) ""; set Cmd(2) "$EUGENEAS $SEQ -s -ph"
set Flag(3) EXECUTION_TRACE3; set Cmd_begin(3) ""; set Cmd(3) "$EUGENEAS $SEQ -s -d"
set Flag(4) EXECUTION_TRACE4; set Cmd_begin(4) ""; set Cmd(4) "$EUGENEAS $SEQ -s -d -E"
set Flag(5) EXECUTION_TRACE5; set Cmd_begin(5) ""; set Cmd(5) "$EUGENEAS $SEQ -s -d -b012"
set Flag(6) EXECUTION_TRACE6; set Cmd_begin(6) "cp $SEQ.user1 $SEQ.user;"; set Cmd(6) "$EUGENEAS $SEQ -U"
set Flag(7) EXECUTION_TRACE7; set Cmd_begin(7) "cp $SEQ.user2 $SEQ.user;"; set Cmd(7) "$EUGENEAS $SEQ -U"
set nbflags 7
#===========================================================================


# read tex doc file
set f [open $FIC_TEX r]
set content [read $f]
close $f


# substitution of flags with execution trace
set new_content ""
for {set i 1} {$i<= $nbflags} {incr i} {
    set FlagPos [string first $Flag($i) $content]
    if { $FlagPos == -1 } {
	puts "ERROR no $Flag($i) set."
	set i $nbflags
    } else {
	set begin_new_content [string range $content 0 [expr $FlagPos - 1]]
	set new_content "$new_content $begin_new_content" 
	set begin_content [expr $FlagPos + [string length $Flag($i)]]
	set content [string range $content $begin_content [string length $content]]

	if {[string length $Cmd_begin($i)] != 0} {
	    eval exec $Cmd_begin($i)
	}

	eval exec $EUGENEAS_DIR$Cmd($i) $Cmd_end

	# write the executed command in the documentation before the result
	set new_content "$new_content >$Cmd($i)\n"
	set f [open $FIC_TMP r]
	set new_content "$new_content [read $f]"
	close $f
    }	
}
set new_content "$new_content $content"

# write new tex doc file
set f [open $FIC_TEX_TMP w]
puts $f $new_content
close $f

# clean directory (beware except image .png)
exec rm EuGeneAS.par  $FIC_TMP

# ask for compilation
exec pdflatex $FIC_TEX_TMP
catch {exec makeindex $CMDFLAGS_INDEX}
exec pdflatex $FIC_TEX_TMP
puts "EuGène documentation generated."
