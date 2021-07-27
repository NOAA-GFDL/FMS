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

# Ryan Mulhall

# Set common test settings.
. ../test-lib.sh

# create and enter directory for in/output
output_dir

#Create file for test.
cat <<_EOF > input.nml
&test_mpp_io_nml
  nx = 360
  ny = 200
  nz = 50
  stackmaxd = 5000000
  layout = 1,1
  io_layout = 1,1
/

&mpp_io_nml
  io_clocks_on = .true.
/
_EOF

test_expect_success "mpp_io functionality with mixed prec reals" '
    mpirun -n 12 ../test_io_R4_R8
'

test_done
