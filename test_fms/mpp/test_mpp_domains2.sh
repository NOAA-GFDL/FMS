#!/bin/sh

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

# This is part of the GFDL FMS package. This is a shell script to
# execute tests in the test_fms/mpp directory.

# Ed Hartnett 11/29/19
# Ryan Mulhall 2/2021

# Set common test settings.
. ../test-lib.sh

# TODO Every test here fails except for 1
SKIP_TESTS="$(basename $0 .sh)."[2,3,4,5,6,7,8,9,10,11,12,13,14,15,16]

# Create base nml for input
sh create_input.sh nml

cat input_base.nml > input.nml # Input.nml is required to run the following tests
test_expect_success "simple functionality" '
    mpirun -n 4 ./test_domains_simple
'
#sed "s/test_nest = .false./test_nest = .true./" $top_srcdir/test_fms/mpp/input_base.nml > input.nml
test_expect_success "update nest domain" '
    mpirun -n 2 ./test_mpp_domains
'
#sed "s/test_subset = .false./test_subset = .true./" $top_srcdir/test_fms/mpp/input_base.nml > input.nml
test_expect_success "subset update" '
    mpirun -n 2 ./test_mpp_domains
'
#sed "s/test_halosize_performance = .false./test_halosize_performance = .true./" $top_srcdir/test_fms/mpp/input_base.nml > input.nml
test_expect_success "halosize performance" '
    mpirun -n 2 ./test_mpp_domains
'
#sed "s/test_edge_update = .false./test_edge_update = .true./" $top_srcdir/test_fms/mpp/input_base.nml > input.nml
test_expect_success "edge update" '
    mpirun -n 2 ./test_mpp_domains
'
#sed "s/test_nonsym_edge = .false./test_nonsym_edge = .true./" $top_srcdir/test_fms/mpp/input_base.nml > input.nml
test_expect_success "nonsym edge" '
    mpirun -n 2 ./test_mpp_domains
'
#sed "s/test_performance = .false./test_performance = .true./" $top_srcdir/test_fms/mpp/input_base.nml > input.nml
#If the system is Darwin or TRAVIS it will be skipped because it fails
test_expect_success "performance" '
    mpirun -n 6 ./test_mpp_domains
'
#sed "s/test_global_sum = .false./test_global_sum = .true./" $top_srcdir/test_fms/mpp/input_base.nml > input.nml
test_expect_success "global sum" '
    mpirun -n 2 ./test_mpp_domains
'
#sed "s/test_cubic_grid_redistribute = .false./test_cubic_grid_redistribute = .true./" $top_srcdir/test_fms/mpp/input_base.nml > input.nml
test_expect_success "cubic grid redistribute" '
    mpirun -n 6 ./test_mpp_domains
'
#sed "s/test_boundary = .false./test_boundary = .true./" $top_srcdir/test_fms/mpp/input_base.nml > input.nml
test_expect_success "boundary" '
    mpirun -n 6 ./test_mpp_domains
'
#sed "s/test_adjoint = .false./test_adjoint = .true./" $top_srcdir/test_fms/mpp/input_base.nml > input.nml
test_expect_success "adjoint" '
    mpirun -n 2 ./test_mpp_domains
'
#sed "s/test_unstruct = .false./test_unstruct = .true./" $top_srcdir/test_fms/mpp/input_base.nml > input.nml
test_expect_success "unstruct" '
    mpirun -n 2 ./test_mpp_domains
'
#sed "s/test_group = .false./test_group = .true./" $top_srcdir/test_fms/mpp/input_base.nml > input.nml
test_expect_success "group" '
    mpirun -n 2 ./test_mpp_domains
'
#sed "s/test_interface = .false./test_interface = .true./" $top_srcdir/test_fms/mpp/input_base.nml > input.nml
test_expect_success "interface" '
    mpirun -n 2 ./test_mpp_domains
'
#sed "s/check_parallel = .false./check_parallel = .true./" $top_srcdir/test_fms/mpp/input_base.nml > input.nml
test_expect_success "check_parallel" '
    mpirun -n 6 ./test_mpp_domains
'
#sed "s/test_get_nbr = .false./test_get_nbr = .true./" $top_srcdir/test_fms/mpp/input_base.nml > input.nml
test_expect_success "get nbr" '
    mpirun -n 2 ./test_mpp_domains
'

test_done
