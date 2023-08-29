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

# 2 bcs * 2 fields = 4 for both dimensions = 8 variables to register
cat <<_EOF > diag_table
test_coupler
1 1 1 0 0 0
#output files
 "coupler_types_bc2",  1, "days", 1, "days", "time"
 "coupler_types_bc1",  1, "days", 1, "days", "time"
#output variables
 "test_coupler_types", "bc1_var2d_1", "bc1_variable_2d_1_min", "coupler_types_bc1", "all", "min", "none", 2
 "test_coupler_types", "bc1_var2d_2", "bc1_variable_2d_2_max", "coupler_types_bc1", "all", "max", "none", 2
 "test_coupler_types", "bc1_var3d_1", "bc1_variable_3d_1", "coupler_types_bc1", "all", "rms", "none", 2
 "test_coupler_types", "bc1_var3d_2", "bc1_variable_3d_2", "coupler_types_bc1", "all", "avg", "none", 2
 "test_coupler_types", "bc2_var2d_1", "bc2_variable_2d_1_min", "coupler_types_bc2", "all", "min", "none", 2
 "test_coupler_types", "bc2_var2d_2", "bc2_variable_2d_2_max", "coupler_types_bc2", "all", "max", "none", 2
 "test_coupler_types", "bc2_var3d_1", "bc2_variable_3d_1", "coupler_types_bc2", "all", "rms", "none", 2
 "test_coupler_types", "bc2_var3d_2", "bc2_variable_3d_2", "coupler_types_bc2", "all", "avg", "none", 2
_EOF

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
