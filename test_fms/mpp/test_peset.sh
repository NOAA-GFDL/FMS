#!/bin/sh

#***********************************************************************
#                   GNU Lesser General Public License
#
# This file is part of the GFDL Flexible Modeling System (FMS).
#
# FMS is free software: you can redistribute it and/or modify it under
# the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or (at
# your option) any later version.
#
# FMS is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
# FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
# for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with FMS.  If not, see <http://www.gnu.org/licenses/>.
#***********************************************************************

# This is part of the GFDL FMS package. This is a shell script to
# execute test_peset in the test_fms/mpp directory.

# Jessica Liptak

# Set common test settings.
. ../test_common.sh

touch input.nml
# Run the test for one processor
echo "Running test_peset with 1 pe"
run_test test_peset 1

# If on a Linux system that uses the command `nproc`, run the test
# with the full number of processors

if [ $(command -v nproc) ]
 # Looks like a linux system
 then
   # Get the number of available CPUs on the system
   nProc=$(nproc)
   if [ ${nProc} -gt 1 ]
     then
       # Run the test with all processors
       echo "Running test_peset with 2 pes"
       run_test test_peset 2
   fi
fi
