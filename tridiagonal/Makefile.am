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

# This is an automake file for the tridiagonal directory of the FMS
# package.

# Ed Hartnett 2/22/19

# Include .h and .mod files.
AM_CPPFLAGS = -I$(top_srcdir)/include -I$(top_srcdir)/tridiagonal/include
AM_FCFLAGS = $(FC_MODINC). $(FC_MODOUT)$(MODDIR)

# Build this uninstalled convenience library.
noinst_LTLIBRARIES = libtridiagonal.la

# The convenience library depends on its source.
libtridiagonal_la_SOURCES = tridiagonal.F90           \
                            include/tridiagonal.inc   \
                            include/tridiagonal_r4.fh \
                            include/tridiagonal_r8.fh

# Mod file depends on its o file, is built and then installed.
tridiagonal.lo: tridiagonal_mod.$(FC_MODEXT)

BUILT_SOURCES = tridiagonal_mod.$(FC_MODEXT)
nodist_include_HEADERS = tridiagonal_mod.$(FC_MODEXT)

include $(top_srcdir)/mkmods.mk
