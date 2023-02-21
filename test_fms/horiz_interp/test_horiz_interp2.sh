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
# Ryan Mulhall 01/23

# Set common test settings.
. ../test-lib.sh

# Create file for test.
cat <<_EOF > input.nml
&test_horiz_interp_nml
  test_conserve = .true.
  ni_src = 360
  nj_src = 180
  ni_dst = 144
  nj_dst = 72
/
_EOF

test_expect_success "conservative method" '
  mpirun -n 2 ./test_horiz_interp_r8
'

test_expect_success "conservative method" '
  mpirun -n 2 ./test_horiz_interp_r4
'

cat <<_EOF > input.nml
&test_horiz_interp_nml
  test_bicubic= .true.
  ni_src = 360
  nj_src = 180
  ni_dst = 144
  nj_dst = 72
/
_EOF

test_expect_success "bicubic method" '
  mpirun -n 2 ./test_horiz_interp
'

cat <<_EOF > input.nml
&test_horiz_interp_nml
  test_bilinear= .true.
  ni_src = 360
  nj_src = 180
  ni_dst = 144
  nj_dst = 72
/
_EOF

test_expect_success "bilinear method" '
  mpirun -n 2 ./test_horiz_interp
'

cat <<_EOF > input.nml
&test_horiz_interp_nml
  test_spherical= .true.
  ni_src = 360
  nj_src = 180
  ni_dst = 360 
  nj_dst = 180
/
_EOF
#TODO
test_expect_success "spherical method" '
  mpirun -n 2 ./test_horiz_interp
'

test_done
