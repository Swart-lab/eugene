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
# Contents: compilation of the User directory
# ------------------------------------------------------------------

AM_CXXFLAGS = $(eugene_CXXFLAGS)
CC = $(CXX)
AM_YFLAGS = -d


bin_PROGRAMS = Sensor.User.so
Sensor_User_so_SOURCES = yacc.y lex.l Sensor.User.cc Sensor.User.h structure.h
Sensor_User_so_LDFLAGS = $(eugene_LDFLAGS)

CLEANFILES = lex.c yacc.h yacc.c *~
