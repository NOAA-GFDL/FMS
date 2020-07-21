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

# MiKyung Lee 06/18/2020, based on Ed Hartnett 11/29/19


# Set common test settings.
. ../test_common.sh


# Run the test with one processor
run_test test_mpp_root_pe 1

# If on a Linux system that uses the command `nproc`,
# Run the test with multiple processors
if [ $(command -v nproc) ]
  # Looks like a linux system
  then
  # Get the number of available CPUs on the system
  nProc=$(nproc)
  if [ ${nProc} -gt 1 ]
  then
    # Run the test with multiple processors
    run_test test_mpp_root_pe 2
  fi
fi
