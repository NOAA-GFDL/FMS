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

test_expect_success "conservative method with real kind=4" '
  mpirun -n 2 ./test_horiz_interp_r8
'

test_expect_success "conservative method with real kind=8" '
  mpirun -n 2 ./test_horiz_interp_r4
'

cat <<_EOF > input.nml
&test_horiz_interp_nml
  test_conserve = .true.
  test_solo = .true.
  ni_src = 360
  nj_src = 180
  ni_dst = 144
  nj_dst = 72
/
_EOF

test_expect_success "conservative method solo wrappers with real kind=4" '
  mpirun -n 2 ./test_horiz_interp_r8
'

test_expect_success "conservative method solo wrappers with real kind=8" '
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

test_expect_success "bicubic method with real kind=4" '
  mpirun -n 2 ./test_horiz_interp_r4
'
test_expect_success "bicubic method with real kind=8" '
  mpirun -n 2 ./test_horiz_interp_r8
'

cat <<_EOF > input.nml
&test_horiz_interp_nml
  test_bicubic= .true.
  test_solo = .true.
  ni_src = 360
  nj_src = 180
  ni_dst = 144
  nj_dst = 72
/
_EOF

test_expect_success "bicubic method solo wrappers with real kind=4" '
  mpirun -n 2 ./test_horiz_interp_r4
'
test_expect_success "bicubic method solo wrappers with real kind=8" '
  mpirun -n 2 ./test_horiz_interp_r8
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

test_expect_success "bilinear method with real kind=4 increasing latitude/longitude" '
  mpirun -n 2 ./test_horiz_interp_r4
'
test_expect_success "bilinear method with real kind=8 increasing latitude/longitude" '
  mpirun -n 2 ./test_horiz_interp_r8
'

cat <<_EOF > input.nml
&test_horiz_interp_nml
  test_bilinear= .true.
  ni_src = 360
  nj_src = 180
  ni_dst = 144
  nj_dst = 72
  decreasing_lat = .true.
/
_EOF

test_expect_success "bilinear method with real kind=4 decreasing latitude/longitude" '
  mpirun -n 2 ./test_horiz_interp_r4
'
test_expect_success "bilinear method with real kind=8 decreasing latitude/longitude" '
  mpirun -n 2 ./test_horiz_interp_r8
'
cat <<_EOF > input.nml
&test_horiz_interp_nml
  test_bilinear= .true.
  test_solo = .true.
  ni_src = 360
  nj_src = 180
  ni_dst = 144
  nj_dst = 72
/
_EOF

test_expect_success "bilinear method solo wrapper with real kind=4 increasing latitude/longitude" '
  mpirun -n 2 ./test_horiz_interp_r4
'
test_expect_success "bilinear method solo wrapper with real kind=8 increasing latitude/longitude" '
  mpirun -n 2 ./test_horiz_interp_r8
'

cat <<_EOF > input.nml
&test_horiz_interp_nml
  test_bilinear= .true.
  test_solo = .true.
  ni_src = 360
  nj_src = 180
  ni_dst = 144
  nj_dst = 72
  decreasing_lat = .true.
/
_EOF

test_expect_success "bilinear method solo wrapper with real kind=4 decreasing latitude/longitude" '
  mpirun -n 2 ./test_horiz_interp_r4
'
test_expect_success "bilinear method solo wrapper with real kind=8 decreasing latitude/longitude" '
  mpirun -n 2 ./test_horiz_interp_r8
'

# the spherical module has a namelist with an option for the search algorithm used
cat <<_EOF > input.nml
&test_horiz_interp_nml
  test_spherical= .true.
  ni_src = 360
  nj_src = 180
  ni_dst = 12
  nj_dst = 6
/

&horiz_interp_sherical_nml
  search_method = "radial search"
/
_EOF

test_expect_success "spherical method (radial search) with real kind=4" '
  mpirun -n 2 ./test_horiz_interp_r4
'
test_expect_success "spherical method (radial search) with real kind=8" '
  mpirun -n 2 ./test_horiz_interp_r8
'

cat <<_EOF > input.nml
&test_horiz_interp_nml
  test_spherical= .true.
  ni_src = 360
  nj_src = 180
  ni_dst = 12
  nj_dst = 6
/

&horiz_interp_sherical_nml
  search_method = "full search"
/
_EOF

test_expect_success "spherical method (full search) with real kind=4" '
  mpirun -n 2 ./test_horiz_interp_r4
'
test_expect_success "spherical method (full search) with real kind=8" '
  mpirun -n 2 ./test_horiz_interp_r8
'

cat <<_EOF > input.nml
&test_horiz_interp_nml
  test_spherical= .true.
  test_solo= .true.
  ni_src = 360
  nj_src = 180
  ni_dst = 12
  nj_dst = 6
/
_EOF

test_expect_success "spherical method solo wrappers with real kind=4" '
  mpirun -n 2 ./test_horiz_interp_r4
'
test_expect_success "spherical method solo wrappers with real kind=8" '
  mpirun -n 2 ./test_horiz_interp_r8
'

cat <<_EOF > input.nml
&test_horiz_interp_nml
  test_assign= .true.
  ni_src = 360
  nj_src = 180
  ni_dst = 12
  nj_dst = 6
/
_EOF

test_expect_success "assignment overloads with real kind=4" '
  mpirun -n 2 ./test_horiz_interp_r4
'
test_expect_success "assignment overloads with real kind=8" '
  mpirun -n 2 ./test_horiz_interp_r8
'
test_expect_success "horiz_interp_read_weights_conserve" 'mpirun -n 4 ./test_horiz_interp_read'

test_expect_success "horiz_interp_fregrid_compatibility" './test_horiz_interp_fregrid_compatibility'

test_done
