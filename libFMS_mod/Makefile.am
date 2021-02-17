#***********************************************************************
#*                   GNU Lesser General Public License
#*
#* This file is part of the GFDL Flexible Modeling System (FMS).
#*
#* FMS is free software: you can redistribute it and/or modify it under
#* the terms of the GNU Lesser General Public License as published by
#* the Free Software Foundation, either version 3 of the License, or (at
#* your option) any later version.
#*
#* FMS is distributed in the hope that it will be useful, but WITHOUT
#* ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
#* FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
#* for more details.
#*
#* You should have received a copy of the GNU Lesser General Public
#* License along with FMS.  If not, see <http://www.gnu.org/licenses/>.
#***********************************************************************

# This is an automake file for the libFMS directory of the FMS
# package. This is the final packaging of the library.

# Ryan Mulhall 2/2021

# This builds the libFMS module file

AM_CPPFLAGS = -I${top_srcdir}/include -I${top_srcdir}/mpp/include
AM_FCFLAGS = $(FC_MODINC). $(FC_MODOUT)$(MODDIR)

noinst_LTLIBRARIES = libFMS.la
libFMS_la_SOURCES = libFMS.F90
libfms_mod.mod: libFMS.$(OBJEXT)
BUILT_SOURCES = libfms_mod.mod
include_HEADERS = libfms_mod.mod
nodist_include_HEADERS = libfms_mod.mod

include $(top_srcdir)/mkmods.mk
