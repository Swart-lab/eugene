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
# File:     Makefile.am
# Contents: compilation of the root directory
# ------------------------------------------------------------------

if LINUX_OS
export sensors_ldflags = -shared
export eugene_CXXFLAGS = -rdynamic
endif

if SOLARIS_OS
export sensors_ldflags = -shared
export eugene_CXXFLAGS = 
endif

if MACOS_OS
export sensors_ldflags = -bundle -flat_namespace -undefined suppres
export eugene_CXXFLAGS = 
endif

EXTRA_DIST = reconf configure  config/* CREATE_LINKS \
	     doc/EuGeneDoc.pdf doc/Sequences/SYNO_ARATH.fasta* \
	     doc/GenerateDoc.tcl doc/Doc.tex doc/Doc_*.tex \
	     doc/ParaOptimization.pdf doc/ParaOptimization.fig \
	     src/getopt.h src/eugene.par \
	     src/PLUGINS/Ara2UTR.mat src/PLUGINS/BLOSUM* src/PLUGINS/WAM/* \
	     src/PLUGINS/swissprot.maxorder2.bin src/PLUGINS/SwP41.noFragm.mininfo1.order2.bin \
	     src/PLUGINS/exon.dist src/PLUGINS/intergenic.dist src/PLUGINS/intron.dist \
	     src/PLUGINS/term.dist src/PLUGINS/init.dist       src/PLUGINS/intr.dist \
             src/PLUGINS/sngl.dist src/PLUGINS/utr.dist \
	     src/WEB/Images/back.jpg src/WEB/Images/euGpred_on.jpg \
             src/WEB/Images/first.jpg src/WEB/Images/last.jpg src/WEB/Images/left.jpg \
	     src/WEB/Javascripts/diap.js \
	     src/WEB/Style/eugene.css 

SUBDIRS = src doc

CLEANFILES = *~

# the files to distribute are automatically found by Automake: 
# all source files are automatically included in a distribution
# the following line avoids to have lex.c and yacc.c 
# in the gzip'd tar file for distribution
dist-hook:
	rm $(distdir)/src/SensorPlugins/User/*.c
