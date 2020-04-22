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

# This is an automake file for the test_fms/drifters directory of the FMS
# package.

# uramirez, Ed Hartnett

# Find the needed mod and .inc files.
AM_CPPFLAGS = -I${top_builddir}/drifters -I${top_builddir}/mpp \
-I${top_srcdir}/include -I${top_builddir}/fms

# Link to the FMS library.
LDADD = ${top_builddir}/libFMS/libFMS.la

# Build these test programs.
check_PROGRAMS = test_cloud_interpolator test_drifters_input		\
test_drifters_comm test_drifters_core test_drifters_io test_quicksort	\
test_drifters

# This is the source code for the tests.
test_cloud_interpolator_SOURCES = test_cloud_interpolator.F90
test_drifters_io_SOURCES = test_drifters_io.F90
test_drifters_input_SOURCES = test_drifters_input.F90
test_drifters_comm_SOURCES = test_drifters_comm.F90
test_drifters_core_SOURCES = test_drifters_core.F90
test_quicksort_SOURCES = test_quicksort.F90
test_drifters_SOURCES = test_drifters.F90

# Run the test programs.
TESTS = test_drifters2.sh

# Include these files in the distribution.
EXTRA_DIST = input_base.nml drifters_inp_test_3d.cdl	\
test_drifters2.sh

# Clean up
CLEANFILES = *out *__genmod.mod *__genmod.f90 input.nml *.nc
