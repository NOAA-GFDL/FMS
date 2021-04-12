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

# TODO 3,4,15

# Create base nml for input
. ./create_input.sh nml

cat input_base.nml > input.nml
test_expect_success "simple functionality" '
    mpirun -n 4 ./test_domains_simple
'
sed "s/test_nest = .false./test_nest = .true./" input_base.nml > input.nml
test_expect_success "update nest domain" '
    mpirun -n 16 ./test_mpp_domains
'
## single face nest needs additional namelist changes
sed "s/tile_coarse =    1,  3,  7/tile_coarse =    1,  1,  2/" input_base.nml > input.nml
sed -i "s/tile_fine   =    7 , 8,  9/tile_fine   =    2 , 3,  4/" input.nml
sed -i "s/istart_coarse =  3,  3,  5/istart_coarse =  4,  3,  5/" input.nml
sed -i "s/icount_coarse = 40,  5,  6/icount_coarse = 12,  5,  6/" input.nml
sed -i "s/jstart_coarse =  3,  3,  6/jstart_coarse =  4,  3,  6/" input.nml
sed -i "s/jcount_coarse = 14,  6,  8/jcount_coarse = 12,  6,  8/" input.nml
sed -i "s/ntiles_nest_all = 9/ntiles_nest_all = 4/" input.nml
sed -i "s/npes_nest_tile = 2, 2, 2, 2, 2, 2, 2, 1, 1/npes_nest_tile = 2, 2, 2, 1/" input.nml
test_expect_success "update single face nest domain" '
    mpirun -n 7 ./test_mpp_domains
'
sed "s/test_subset = .false./test_subset = .true./" input_base.nml > input.nml
test_expect_success "subset update" '
    mpirun -n 25 ./test_mpp_domains
'
sed "s/test_halosize_performance = .false./test_halosize_performance = .true./" input_base.nml > input.nml
test_expect_success "halosize performance" '
    mpirun -n 2 ./test_mpp_domains
'
sed "s/test_edge_update = .false./test_edge_update = .true./" input_base.nml > input.nml
test_expect_success "edge update" '
    mpirun -n 2 ./test_mpp_domains
'
sed "s/test_nonsym_edge = .false./test_nonsym_edge = .true./" input_base.nml > input.nml
test_expect_success "nonsym edge" '
    mpirun -n 2 ./test_mpp_domains
'
sed "s/test_performance = .false./test_performance = .true./" input_base.nml > input.nml
#If the system is Darwin or TRAVIS it will be skipped because it fails
test_expect_success "performance" '
    mpirun -n 6 ./test_mpp_domains
'
sed "s/test_global_sum = .false./test_global_sum = .true./" input_base.nml > input.nml
test_expect_success "global sum" '
    mpirun -n 2 ./test_mpp_domains
'
sed "s/test_cubic_grid_redistribute = .false./test_cubic_grid_redistribute = .true./" input_base.nml > input.nml
test_expect_success "cubic grid redistribute" '
    mpirun -n 6 ./test_mpp_domains
'
sed "s/test_boundary = .false./test_boundary = .true./" input_base.nml > input.nml
test_expect_success "boundary" '
    mpirun -n 6 ./test_mpp_domains
'
sed "s/test_adjoint = .false./test_adjoint = .true./" input_base.nml > input.nml
test_expect_success "adjoint" '
    mpirun -n 2 ./test_mpp_domains
'
sed "s/test_unstruct = .false./test_unstruct = .true./" input_base.nml > input.nml
test_expect_success "unstruct" '
    mpirun -n 2 ./test_mpp_domains
'
sed "s/test_group = .false./test_group = .true./" input_base.nml > input.nml
test_expect_success "group" '
    mpirun -n 2 ./test_mpp_domains
'
sed "s/test_interface = .false./test_interface = .true./" input_base.nml > input.nml
test_expect_success "interface" '
    mpirun -n 2 ./test_mpp_domains
'
# TODO this wasn't run at all previously
#sed "s/check_parallel = .false./check_parallel = .true./" input_base.nml > input.nml
#test_expect_success "check_parallel" '
#    mpirun -n 2 ./test_mpp_domains
#'
sed "s/test_get_nbr = .false./test_get_nbr = .true./" input_base.nml > input.nml
test_expect_success "get nbr" '
    mpirun -n 8 ./test_mpp_domains
'

test_done
