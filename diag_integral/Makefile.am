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

# This is an automake file for the diag_integral directory of the FMS
# package.

# Ed Hartnett 2/22/19

# Include .h and .mod files.
AM_CPPFLAGS = -I$(top_srcdir)/include -I$(top_srcdir)/diag_integral/include
AM_FCFLAGS = $(FC_MODINC). $(FC_MODOUT)$(MODDIR)

# Build this uninstalled convenience library.
noinst_LTLIBRARIES = libdiag_integral.la

# The convenience library depends on its source.
libdiag_integral_la_SOURCES = diag_integral.F90\
															include/diag_integral.inc\
															include/diag_integral_r4.fh\
															include/diag_integral_r8.fh

diag_integral_mod.$(FC_MODEXT): include/diag_integral.inc\
																include/diag_integral_r4.fh\
																include/diag_integral_r8.fh


nodist_include_HEADERS = diag_integral_mod.$(FC_MODEXT)
BUILT_SOURCES = diag_integral_mod.$(FC_MODEXT)


include $(top_srcdir)/mkmods.mk
