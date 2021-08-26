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
# execute tests in the test_fms/horiz_interp directory.

# Ed Hartnett 11/29/19

# Set common test settings.
. ../test-lib.sh

# Create file for test.
cat <<_EOF > input.nml
&test_horiz_interp_nml
  ni_src = 360
  nj_src = 180
  ni_dst = 144
  nj_dst = 72
/
_EOF

test_expect_success "Horiz_interp test" '
  mpirun -n 2 ./test_horiz_interp
'

test_done
