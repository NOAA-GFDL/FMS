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

# This is an automake file for the tracer_manager directory of the FMS
# package.

# Ed Hartnett 2/22/19

# Include .h and .mod files.
AM_CPPFLAGS = -I$(top_srcdir)/include -I$(top_srcdir)/tracer_manager/include
AM_FCFLAGS = $(FC_MODINC). $(FC_MODOUT)$(MODDIR)

# Build this uninstalled convenience library.
noinst_LTLIBRARIES = libtracer_manager.la

# The convenience library depends on its source.
libtracer_manager_la_SOURCES = tracer_manager.F90 \
															 include/tracer_manager.inc\
															 include/tracer_manager_r4.fh\
															 include/tracer_manager_r8.fh


BUILT_SOURCES = tracer_manager_mod.$(FC_MODEXT) include/tracer_manager.inc include/tracer_manager_r4.fh include/tracer_manager_r8.fh
nodist_include_HEADERS = tracer_manager_mod.$(FC_MODEXT)

include $(top_srcdir)/mkmods.mk
