# ------------------------------------------------------------------
# Copyright (C) 2004 INRA <eugene-help@mulcyber.toulouse.inra.fr>
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
# Contents: compilation of the SignalWAM directory
# ------------------------------------------------------------------

AM_CXXFLAGS = $(eugene_cxxflags)
AM_CFLAGS = 


noinst_PROGRAMS = Sensor.SignalWAM.so
Sensor_SignalWAM_so_SOURCES = Sensor.SignalWAM.cc Sensor.SignalWAM.h 
Sensor_SignalWAM_so_LDFLAGS = $(sensors_ldflags) 


CLEANFILES = *~
