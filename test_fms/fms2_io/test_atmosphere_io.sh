#!/bin/sh

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
