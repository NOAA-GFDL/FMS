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

# Author: Ed Hartnett 6/10/20
#
# Set common test settings.
. ../test-lib.sh

# Create and enter output directory
output_dir

# make an input.nml for mpp_init to read
touch input.nml

# run the tests
test_expect_success "Test the filename_appendix functionality" '
  mpirun -n 1 ../test_file_appendix
'
test_expect_success "Simple IO test" '
  mpirun -n 6 ../test_io_simple
'
test_expect_success "Test the get_mosaic_tile_grid functionality" '
  mpirun -n 6 ../test_get_mosaic_tile_grid
'
test_expect_success "Test the get_valid is_valid functionality single PE" '
  mpirun -n 1 ../test_get_is_valid
'
test_expect_success "Test the get_valid is_valid functionality multiple PE" '
  mpirun -n 1 ../test_get_is_valid
'
test_expect_success "Test the unlimited compressed axis functionality" '
  mpirun -n 6 ../test_unlimit_compressed
'

test_done
