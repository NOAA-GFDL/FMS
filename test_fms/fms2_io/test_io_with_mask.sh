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

# Author: Uriel Ramirez 07/07/20
#
# Set common test settings.
. ../test_common.sh

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
printf "\n1\n2,3\n1,1" | cat > the_mask

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
run_test test_io_with_mask 5
