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
# execute tests in the test_fms/mosaic directory.

# Ed Hartnett 11/29/19

# Set common test settings.
. ../test-lib.sh

# Copy files for test.
touch input.nml
rm -rf INPUT
mkdir INPUT

if [ ! $parser_skip ]; then
  SKIP_TESTS='test_mosaic2.[1-4]'
fi

# The tests are skipped if FMS is compiled in r4 via ./configure --enable-mixedmode
# because answers differ when FMS is compiled in r4.
test_expect_success "test mosaic2 r4" 'mpirun -n 1 ./test_mosaic2_r4'
test_expect_success "test grid2   r4" 'mpirun -n 1 ./test_grid2_r4'
test_expect_success "test mosaic2 r8" 'mpirun -n 1 ./test_mosaic2_r8'
test_expect_success "test grid2   r8" 'mpirun -n 1 ./test_grid2_r8'

rm -rf INPUT
test_done
