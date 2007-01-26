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



export sensors_ldflags = -shared 
export eugene_cxxflags = -rdynamic -fPIC
export eugene_cflags = -fPIC

if SOLARIS_OS
export sensors_ldflags = -shared -mimpure-text
export eugene_cxxflags = 
export eugene_cflags = 
endif

if MACOS_OS
export sensors_ldflags = -bundle -flat_namespace -undefined suppress
export eugene_cxxflags = 
export eugene_cflags = 
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

# note that $(pkgdatadir) = $(datadir)/$(PACKAGE)
install-data-hook:
	$(INSTALL) -d $(DESTDIR)/$(pkgdatadir)
	$(INSTALL) -d $(DESTDIR)/$(pkgdatadir)/doc 
	$(INSTALL) -d $(DESTDIR)/$(pkgdatadir)/doc/Sequences 
	$(INSTALL) -d $(DESTDIR)/$(pkgdatadir)/web
	$(INSTALL) -d $(DESTDIR)/$(pkgdatadir)/web/Style
	$(INSTALL) -d $(DESTDIR)/$(pkgdatadir)/web/Javascripts
	$(INSTALL) -d $(DESTDIR)/$(pkgdatadir)/web/Images
	$(INSTALL) -d $(DESTDIR)/$(pkgdatadir)/plugins
	$(INSTALL) -d $(DESTDIR)/$(pkgdatadir)/cfg
	$(INSTALL) -d $(DESTDIR)/$(pkgdatadir)/models
	$(INSTALL) -d $(DESTDIR)/$(pkgdatadir)/models/WAM
	$(INSTALL) -d $(DESTDIR)/$(pkgdatadir)/Procedures
	$(INSTALL) -d $(DESTDIR)/$(pkgdatadir)/Procedures/Eval
	$(INSTALL) -d $(DESTDIR)/$(pkgdatadir)/Procedures/Get
	$(INSTALL) -m 644 $(srcdir)/doc/EuGeneDoc.pdf $(DESTDIR)/$(pkgdatadir)/doc
	$(INSTALL) -m 644 $(srcdir)/doc/Sequences/SYNO_ARATH.fasta* $(DESTDIR)/$(pkgdatadir)/doc/Sequences
	$(INSTALL) -m 644 $(srcdir)/models/*.dist $(DESTDIR)/$(pkgdatadir)/models
	$(INSTALL) -m 644 $(srcdir)/models/*.mat $(DESTDIR)/$(pkgdatadir)/models
	$(INSTALL) -m 644 $(srcdir)/models/*.bin $(DESTDIR)/$(pkgdatadir)/models
	$(INSTALL) -m 644 $(srcdir)/models/BLOSUM* $(DESTDIR)/$(pkgdatadir)/models
	$(INSTALL) -m 644 $(srcdir)/models/WAM/WAM* $(DESTDIR)/$(pkgdatadir)/models/WAM
	$(INSTALL) -m 644 $(srcdir)/web/Style/eugene.css $(DESTDIR)/$(pkgdatadir)/web/Style
	$(INSTALL) -m 644 $(srcdir)/web/Javascripts/diap.js $(DESTDIR)/$(pkgdatadir)/web/Javascripts
	$(INSTALL) -m 644 $(srcdir)/web/Images/*jpg $(DESTDIR)/$(pkgdatadir)/web/Images
	$(INSTALL) -m 644 $(srcdir)/cfg/*.par $(DESTDIR)/$(pkgdatadir)/cfg
	$(INSTALL) src/SensorPlugins/*/*.so  $(DESTDIR)/$(pkgdatadir)/plugins
	$(INSTALL) $(srcdir)/Procedures/Eval/egn_* $(DESTDIR)/$(pkgdatadir)/Procedures/Eval
	$(INSTALL) $(srcdir)/Procedures/Get/egn_* $(DESTDIR)/$(pkgdatadir)/Procedures/Get
