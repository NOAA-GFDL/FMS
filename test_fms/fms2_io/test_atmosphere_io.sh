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
# execute tests in the test_fms/fms2_io directory.

# Set common test settings.
. ../test_common.sh

[ -d atmosphere-output ] && rm -r atmosphere-output

# run test 1 - standard test
mkdir atmosphere-output
cd atmosphere-output
touch input.nml
run_test ../test_atmosphere_io 6
cd .. && rm -r atmosphere-output

# run test 2 - test for bad checksum (should fail)
mkdir atmosphere-output
cd atmosphere-output
printf "&test_atmosphere_io_nml\n bad_checksum=.true.\n /" | cat > input.nml
if run_test ../test_atmosphere_io 6; then
  echo " test was supposed to fail and didn't "
  exit -99
fi
cd .. && rm -r atmosphere-output

# run test 3 - test for ignoring a bad checksum
mkdir atmosphere-output
cd atmosphere-output
printf "&test_atmosphere_io_nml\n bad_checksum=.true.\n ignore_checksum=.true.\n /" | cat > input.nml
run_test ../test_atmosphere_io 6
cd .. && rm -r atmosphere-output

