#!/bin/sh

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

# This is part of the GFDL FMS package. This is a shell script to
# execute tests in the test_fms/mpp directory.

# Ed Hartnett 11/29/19

# Set common test settings.
. ../test_common.sh

skip_test="no"

# Copy file for test.
cp $top_srcdir/test_fms/mpp_io/input_base.nml input.nml


# Get the number of available CPUs on the system
if [ $(command -v nproc) ]
then
   # Looks like a linux system
   nProc=$(nproc)
elif [ $(command -v sysctl) ]
then
   # Looks like a Mac OS X system
   nProc=$(sysctl -n hw.physicalcpu)
else
   nProc=-1
fi

if [ $nProc -lt 0 ]
then
   # Couldn't get the number of CPUs, skip the test.
   skip_test="skip"
elif [ $nProc -lt 12 ]
then
   # Need to oversubscribe the MPI
   run_test test_mpp_io 12 $skip_test "true"
fi

run_test test_mpp_io 12 $skip_test
