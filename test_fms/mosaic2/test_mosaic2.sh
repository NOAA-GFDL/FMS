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
