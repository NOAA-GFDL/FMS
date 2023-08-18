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
# execute tests in the test_fms/data_override directory.

# Ed Hartnett 11/26/19
# Uriel Ramirez 07/22/20

# Set common test settings.
. ../test-lib.sh

# Run the ongrid test case with 2 halos in x and y
touch input.nml

mkdir RESTART

test_expect_success "coupler register restart 2D(r4_kind)" '
  mpirun -n 1 ./test_coupler_2d_r4
'
test_expect_success "coupler register restart 2D(r8_kind)" '
  mpirun -n 1 ./test_coupler_2d_r8
'

test_expect_success "coupler register restart 3D (r4_kind)" '
  mpirun -n 1 ./test_coupler_3d_r4
'

test_expect_success "coupler register restart 3D (r8_kind)" '
  mpirun -n 1 ./test_coupler_3d_r8
'

test_expect_success "coupler type operation interfaces (r4_kind)" '
  mpirun -n 4 ./test_coupler_types_r4
'

test_expect_success "coupler type operation interfaces (r8_kind)" '
  mpirun -n 4 ./test_coupler_types_r8
'

test_expect_success "test atmos_ocean_fluxes" '
  mpirun -n 1 ./test_atmos_ocean_fluxes
'

rm -rf RESTART

test_done
