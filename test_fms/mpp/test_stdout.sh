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
# execute tests in the test_fms/mpp directory.

# Colin Gladue 06/12/20

# Set common test settings.
. ../test_common.sh

# Run test with one processor
echo "Running stdout unit test with 1 proc..."; echo
run_test test_stdout 1
echo; echo "stdout test passed with 1 proc."; echo

# If on a Linux system that uses the command `nproc`, run the test
# with multiple processors

if [ $(command -v nproc) ]
  # Looks like a linux system
  then
  # Get the number of available CPUs on the system
  nProc=$(nproc)
  if [ ${nProc} -gt 1 ]
  then
    # Run the test with multiple processors
    echo "Running stdout unit test with multiple procs..."; echo
    err=0
    run_test test_stdout 2 || err=1
    if [ $err -eq 1 ]; then
      echo; echo "Failed with multiple procs"
      exit 1
    fi
    echo; echo "stdout test passed with multiple procs."
  fi
fi

