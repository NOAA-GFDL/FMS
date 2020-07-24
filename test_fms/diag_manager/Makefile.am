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

# This is an automake file for the test_fms/diag_manager directory of the FMS
# package.
# uramirez, Ed Hartnett

# Find the needed mod and .inc files.
AM_CPPFLAGS = -I${top_srcdir}/include \
  -I${top_builddir}/mpp	\
  -I${top_builddir}/fms \
  -I${top_builddir}/axis_utils \
  -I${top_builddir}/diag_manager \
  -I${top_builddir}/time_manager \
  -I${top_builddir}/constants \
  -I${top_builddir}/platform

# Link to the FMS library.
LDADD = ${top_builddir}/libFMS/libFMS.la

# Build this test program.
check_PROGRAMS = test_diag_manager

# This is the source code for the test.
test_diag_manager_SOURCES = test_diag_manager.F90

# Run the test.
TESTS = test_diag_manager2.sh

# Copy over other needed files to the srcdir
EXTRA_DIST = input.nml_base diagTables test_diag_manager2.sh

CLEANFILES = input.nml *.nc *.out diag_table
