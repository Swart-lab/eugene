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
# Contents: compilation of the 0_SensorTk directory
# Comments: The library is not useful, but allows to automatically define source files for distribution.
# ------------------------------------------------------------------

CXXFLAGS = $(eugene_CXXFLAGS)

SUBDIRS = GetData

noinst_LIBRARIES = libSensorTk.a
libSensorTk_a_SOURCES = BStrArray.cc BStrArray.h EndianConv.h WAM.cc WAM.h markov.cc markov.h



CLEANFILES = *~
