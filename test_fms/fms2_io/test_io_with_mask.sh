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
# execute tests in the test_fms/fms2_io directory.

# Author: Uriel Ramirez 07/07/20
#
# Set common test settings.
. ../test-lib.sh

# Create and enter output directory
output_dir

# make an input.nml for mpp_init to read
touch input.nml

# make a masktable
# This mask table masks out 1 rank (1st line), for a layout of 2,3 (2nd line)
# 1,1 is the section that gets masked out
# . ----- . ----- .
# | (1,1) | (1,2) |
# . ----- . ----- .
# | (2,1) | (2,2) |
# . ----- . ----- .
# | (3,1) | (3,2) |
# . ----- . ----- .

printf "1\n2,3\n1,1" | cat > the_mask

# For example, if you have a grid that is 60 by 60 and a layout of 2,3
# You are going to need 6 ranks:
   # rank 1 is going to handle: x: 1-30 y: 1-20
   # rank 2 is going to handle: x: 31-60 y: 1-20
   # rank 3 is going to handle: x: 1-30 y: 21-40
   # rank 4 is going to handle: x: 31-60 y: 21-40
   # rank 5 is going to handle: x: 1-30 y: 41-60
   # rank 6 is going to handle: x: 31-60 y: 41-60
# But with the mask table above, rank 0 is going to be masked out so you know need
# 5 ranks and nothing is going to be done for x: 1-30 y: 1-20

# run the tests
test_expect_success "Test FMS2 IO using a mask" '
  mpirun -n 5 ../test_io_with_mask
'

test_done
