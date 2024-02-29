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
. ../test-lib.sh

# Create and enter output directory
output_dir

touch input.nml

# run the tests
test_expect_success "Test atmosphere IO" '
  mpirun -n 6 ../test_atmosphere_io
'

# run test 2 - test for bad checksum (should fail)
rm *.nc
printf "&test_atmosphere_io_nml\n bad_checksum=.true.\n /" | cat > input.nml
test_expect_failure "bad checksum failure" '
  mpirun -n 6 ../test_atmosphere_io
'

# run test 3 - test for ignoring a bad checksum
rm *.nc
printf "&test_atmosphere_io_nml\n bad_checksum=.true.\n ignore_checksum=.true.\n /" | cat > input.nml
test_expect_success "ignore bad checksum" '
  mpirun -n 6 ../test_atmosphere_io
'

test_done
