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
# Contents: compilation of the src directory
# ------------------------------------------------------------------

AUTOMAKE_OPTIONS = dejagnu


SUBDIRS = Parametrization GDIF . SensorPlugins

CXXFLAGS = $(eugene_CXXFLAGS)

LIBS = `if [ "$$MACHTYPE" = "sparc" ] ; then \
	echo "../lib/gd-2.0.22/libgd.a ../lib/libpng-1.2.5/libpng.a ../lib/zlib-1.2.1/libz.a"; \
else \
	echo "-lgd -lpng -lz"; \
fi`

bin_PROGRAMS = eugene
eugene_SOURCES = EuGene.cc Output.cc Output.h Hits.cc Hits.h DNASeq.cc DNASeq.h  System.cc System.h Prediction.cc Prediction.h  BackP.cc BackP.h Sensor.cc Sensor.h Dll.cc Dll.h MSensor.cc MSensor.h Param.cc Param.h PenaltyDist.cc PenaltyDist.h SensorIF.h Const.h
eugene_LDADD = Parametrization/libParametrization.a  GDIF/libGDIF.a -ldl -lnsl


TEST = 
RUNTESTDEFAULTFLAGS =  --tool TestSuite --srcdir ../Test --objdir . --outdir ../Test/TestTrace $(TEST)

CLEANFILES = *~
