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
export eugene_cxxflags = -rdynamic
endif

if SOLARIS_OS
export sensors_ldflags = -shared -mimpure-text
export eugene_cxxflags = 
endif

if MACOS_OS
export sensors_ldflags = -bundle -flat_namespace -undefined suppres
export eugene_cxxflags = 
endif


SUBDIRS = src plugins doc 

CLEANFILES = *~


############## Tests
AUTOMAKE_OPTIONS = dejagnu

export host_alias = @host@
TEST = 
RUNTESTDEFAULTFLAGS =  --tool TestSuite --srcdir Test --objdir . --outdir Test/TestTrace $(TEST)


############## Distribute
EXTRA_DIST = reconf configure $(srcdir)/config/config.* \
	     $(srcdir)/config/install-sh \
             $(srcdir)/config/m* \
	     $(srcdir)/doc/EuGeneDoc.pdf \
	     $(srcdir)/doc/Sequences/SYNO_ARATH.fasta* \
	     $(srcdir)/doc/GenerateDoc.tcl \
	     $(srcdir)/doc/Doc.tex \
	     $(srcdir)/doc/Doc_*.tex \
	     $(srcdir)/doc/ParaOptimization.pdf \
	     $(srcdir)/doc/ParaOptimization.fig \
	     $(srcdir)/src/getopt.h \
	     $(srcdir)/cfg/eugene.par \
	     $(srcdir)/models/*.dist \
	     $(srcdir)/models/*.mat \
	     $(srcdir)/models/*.bin \
	     $(srcdir)/models/BLOSUM* \
	     $(srcdir)/models/WAM/WAM* \
             $(srcdir)/web/Style/eugene.css \
             $(srcdir)/web/Javascripts/diap.js \
             $(srcdir)/web/Images/*.jpg \
             $(srcdir)/Procedures/Eval/egn_evalpred.pl \
             $(srcdir)/Procedures/Get/egn_*

# the files to distribute are automatically found by Automake: 
# all source files are automatically included in a distribution
# the following line avoids to have lex.c and yacc.c 
# in the gzip'd tar file for distribution
dist-hook:
	if [ -f $(distdir)/src/SensorPlugins/User/lex.c ]; then rm $(distdir)/src/SensorPlugins/User/lex.c; fi
	if [ -f $(distdir)/src/SensorPlugins/User/yacc.c ]; then rm $(distdir)/src/SensorPlugins/User/yacc.c; fi




############## Install / Uninstall
uninstall-local:
	rm -r $(datadir)/$(PACKAGE)

install-data-hook:
	# note that $(pkgdatadir) = $(datadir)/$(PACKAGE)
	if ! [ -d $(datadir) ]; then mkdir $(datadir); fi
	mkdir $(pkgdatadir)
	mkdir $(pkgdatadir)/doc 
	mkdir $(pkgdatadir)/doc/Sequences 
	mkdir $(pkgdatadir)/web
	mkdir $(pkgdatadir)/web/Style
	mkdir $(pkgdatadir)/web/Javascripts
	mkdir $(pkgdatadir)/web/Images
	mkdir $(pkgdatadir)/plugins
	mkdir $(pkgdatadir)/cfg
	mkdir $(pkgdatadir)/models
	mkdir $(pkgdatadir)/models/WAM
	mkdir $(pkgdatadir)/Procedures
	mkdir $(pkgdatadir)/Procedures/Eval
	mkdir $(pkgdatadir)/Procedures/Get
	cp $(srcdir)/doc/EuGeneDoc.pdf $(pkgdatadir)/doc
	cp $(srcdir)/doc/Sequences/SYNO_ARATH.fasta* $(pkgdatadir)/doc/Sequences
	cp $(srcdir)/models/*.dist $(pkgdatadir)/models
	cp $(srcdir)/models/*.mat $(pkgdatadir)/models
	cp $(srcdir)/models/*.bin $(pkgdatadir)/models
	cp $(srcdir)/models/BLOSUM* $(pkgdatadir)/models
	cp $(srcdir)/models/WAM/WAM* $(pkgdatadir)/models/WAM
	cp $(srcdir)/web/Style/eugene.css $(pkgdatadir)/web/Style
	cp $(srcdir)/web/Javascripts/diap.js $(pkgdatadir)/web/Javascripts
	cp $(srcdir)/web/Images/*.jpg $(pkgdatadir)/web/Images
	cp $(srcdir)/cfg/eugene.par $(pkgdatadir)/cfg
	cp src/SensorPlugins/*/*.so $(pkgdatadir)/plugins
	cp $(srcdir)/Procedures/Eval/egn_evalpred.pl $(pkgdatadir)/Procedures/Eval
	cp $(srcdir)/Procedures/Get/egn_* $(pkgdatadir)/Procedures/Get
