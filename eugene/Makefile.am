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
export sensors_ldflags = -bundle -flat_namespace -undefined suppres
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
	install -d $(DESTDIR)/$(pkgdatadir)
	install -d $(DESTDIR)/$(pkgdatadir)/doc 
	install -d $(DESTDIR)/$(pkgdatadir)/doc/Sequences 
	install -d $(DESTDIR)/$(pkgdatadir)/web
	install -d $(DESTDIR)/$(pkgdatadir)/web/Style
	install -d $(DESTDIR)/$(pkgdatadir)/web/Javascripts
	install -d $(DESTDIR)/$(pkgdatadir)/web/Images
	install -d $(DESTDIR)/$(pkgdatadir)/plugins
	install -d $(DESTDIR)/$(pkgdatadir)/cfg
	install -d $(DESTDIR)/$(pkgdatadir)/models
	install -d $(DESTDIR)/$(pkgdatadir)/models/WAM
	install -d $(DESTDIR)/$(pkgdatadir)/Procedures
	install -d $(DESTDIR)/$(pkgdatadir)/Procedures/Eval
	install -d $(DESTDIR)/$(pkgdatadir)/Procedures/Get
	install -m a-x $(srcdir)/doc/EuGeneDoc.pdf $(DESTDIR)/$(pkgdatadir)/doc
	install -m a-x $(srcdir)/doc/Sequences/SYNO_ARATH.fasta* $(DESTDIR)/$(pkgdatadir)/doc/Sequences
	install -m a-x $(srcdir)/models/*.dist $(DESTDIR)/$(pkgdatadir)/models
	install -m a-x $(srcdir)/models/*.mat $(DESTDIR)/$(pkgdatadir)/models
	install -m a-x $(srcdir)/models/*.bin $(DESTDIR)/$(pkgdatadir)/models
	install -m a-x $(srcdir)/models/BLOSUM* $(DESTDIR)/$(pkgdatadir)/models
	install -m a-x $(srcdir)/models/WAM/WAM* $(DESTDIR)/$(pkgdatadir)/models/WAM
	install -m a-x $(srcdir)/web/Style/eugene.css $(DESTDIR)/$(pkgdatadir)/web/Style
	install -m a-x $(srcdir)/web/Javascripts/diap.js $(DESTDIR)/$(pkgdatadir)/web/Javascripts
	install -m a-x $(srcdir)/web/Images/*jpg $(DESTDIR)/$(pkgdatadir)/web/Images
	install -m a-x $(srcdir)/cfg/*.par $(DESTDIR)/$(pkgdatadir)/cfg
	install src/SensorPlugins/*/*.so  $(DESTDIR)/$(pkgdatadir)/plugins
	install $(srcdir)/Procedures/Eval/egn_* $(DESTDIR)/$(pkgdatadir)/Procedures/Eval
	install $(srcdir)/Procedures/Get/egn_* $(DESTDIR)/$(pkgdatadir)/Procedures/Get
