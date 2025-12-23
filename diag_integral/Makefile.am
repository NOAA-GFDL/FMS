#***********************************************************************
#*                             Apache License 2.0
#*
#* This file is part of the GFDL Flexible Modeling System (FMS).
#*
#* Licensed under the Apache License, Version 2.0 (the "License");
#* you may not use this file except in compliance with the License.
#* You may obtain a copy of the License at
#*
#*     http://www.apache.org/licenses/LICENSE-2.0
#*
#* FMS is distributed in the hope that it will be useful, but WITHOUT
#* WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied;
#* without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
#* PARTICULAR PURPOSE. See the License for the specific language
#* governing permissions and limitations under the License.
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
