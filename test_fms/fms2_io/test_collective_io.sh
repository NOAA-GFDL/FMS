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

#
# Set common test settings.
. ../test-lib.sh

if [ ! -z $parallel_skip ]; then
  SKIP_TESTS="test_collective_io.[1-3]"
fi

# Create and enter output directory
output_dir

touch input.nml

test_expect_success "Test NetCDF-4 MPI-IO parallel writes" '
  mpirun -n 6 ../test_parallel_writes
'

test_expect_success "Test NetCDF-4 MPI-IO collective reads" '
  mpirun -n 6 ../test_collective_io
'

rm -rf *.nc*
cat <<_EOF > input.nml
&test_collective_io_nml
  nc_format = "64bit"
/
_EOF

test_expect_failure "Attempt to open a 64-bit NetCDF file for MPI-IO collective reads" '
  mpirun -n 6 ../test_collective_io
'

test_done
