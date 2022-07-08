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

# Authors: Raymond Menzel
# Jessica Liptak
#
# Set common test settings.
. ../test-lib.sh

# Create and enter output directory
output_dir

# use smaller arrays if system stack size is limited
if [ $STACK_LIMITED ]; then
  cat <<_EOF > input.nml
&test_fms2_io_nml
  nx = 32
  ny = 32
  nz = 10
/
_EOF
fi
touch input.nml

# run the tests
test_expect_success "FMS2 IO Test" '
  mpirun -n 6 ../test_fms2_io
'

test_done
