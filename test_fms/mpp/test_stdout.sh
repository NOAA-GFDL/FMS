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


# If more than one processor available:
echo "Running stdout unit test with multiple procs..."; echo
if [ $(command -v nproc) ]
  then
    nProc=$(nproc)
    run_test test_stdout ${nProc}
fi
echo; echo "stdout test passed with multiple procs."
