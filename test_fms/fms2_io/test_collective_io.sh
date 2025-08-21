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

#
# Set common test settings.
. ../test-lib.sh

if [ ! -z $parallel_skip ]; then
  SKIP_TESTS="test_collective_io.[1-2]"
fi

# Create and enter output directory
output_dir

touch input.nml

test_expect_success "Test NetCDF-4 parallel writes" '
  mpirun -n 6 ../test_parallel_writes
'

test_expect_success "Test NetCDF-4 collective reads" '
  mpirun -n 6 ../test_collective_io
'

rm -rf *.nc*
# The code should still run if not using netcdf4 files, it just won't use collective io
cat <<_EOF > input.nml
&test_collective_io_nml
  nc_format = "64bit"
/
_EOF

test_expect_success "Test fallback to non-collective reads" '
  mpirun -n 6 ../test_collective_io
'

test_done
