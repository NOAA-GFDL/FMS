#!/bin/sh
#
# author @underwoo
#
#***********************************************************************
#                   GNU Lesser General Public License
#
# This file is part of the GFDL Flexible Modeling System (FMS).
#
# FMS is free software: you can redistribute it and/or modify it under
# the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or (at
# your option) any later version.
#
# FMS is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
# FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
# for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with FMS.  If not, see <http://www.gnu.org/licenses/>.
#***********************************************************************
#
# Script to test the mpp_memutils_mod Fortran module code.

# Set common test settings
. ../test_common.sh

# All tests use blank input.nml file.
touch input.nml

echo "1: Test begin/end routines of mpp_memutils_mod"
run_test ./test_mpp_memutils_begin_end 1

echo "2: Test mpp_print_memuse_stats"
run_test ./test_mpp_print_memuse_stats_stderr 1

echo "3: Test mpp_print_memuse_stats to file (stdout)"
run_test ./test_mpp_print_memuse_stats_file 1

echo "4: Test failure caught if mpp_memuse_begin called multiple times"
run_test ./test_mpp_memutils_begin_2x 1 || echo "Ok"

echo "5: Test failure caught if mpp_memuse_end called before mpp_memuse_begin"
run_test ./test_mpp_memutils_end_before_begin 1 || echo "Ok"
