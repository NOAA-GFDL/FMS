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

# This is an automake file for the global directory of the FMS
# package.
#
# This directory only holds the libFMS.F90 file, which gets compiled into
# a module just named 'fms'. It is meant to access all public routines/types/variables from
# the other modules in fms.

# Ryan Mulhall 1/19/24

AM_CPPFLAGS = -I${top_srcdir}/include
AM_FCFLAGS = $(FC_MODINC). $(FC_MODOUT)$(MODDIR)

noinst_LTLIBRARIES = libfms_global.la

libfms_global_la_SOURCES = fms_global.F90

BUILT_SOURCES = fms_global_mod.$(FC_MODEXT)
nodist_include_HEADERS = fms_global_mod.$(FC_MODEXT)

include $(top_srcdir)/mkmods.mk
