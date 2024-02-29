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

# TODO edge update, fails on non-blocking with gnu
#SKIP_TESTS="$SKIP_TESTS $(basename $0 .sh).6"

# create and enter directory for in/output
output_dir

# Create input_base.nml for test input
cat <<_EOF > input_base.nml
&test_mpp_domains_nml
nx=64
ny=64
nz=10
stackmax=10000000
debug=.false.
mpes = 3
check_parallel = .false.
whalo = 2
ehalo = 2
shalo = 2
nhalo = 2
x_cyclic_offset = 3
y_cyclic_offset = -4
warn_level = "fatal"
wide_halo_x = 0
wide_halo_y = 0
nx_cubic = 20
ny_cubic = 20
test_performance = .false.
test_interface = .false.
num_fields = 4
do_sleep = .false.
num_iter = 1
mix_2D_3D = .false.
test_get_nbr = .false.
test_edge_update = .false.
test_cubic_grid_redistribute = .false.
ensemble_size = 1
layout_cubic = 0,0
layout_ensemble = 0,0
nthreads = 1
test_boundary = .false.
layout_tripolar = 0,0
test_group = .false.
test_global_sum = .false.
test_unstruct = .false.
test_nonsym_edge = .false.
test_halosize_performance = .false.
test_adjoint = .false.
wide_halo = .false.
/
_EOF

cat input_base.nml > input.nml
test_expect_success "simple functionality" '
    mpirun -n 4 ../test_domains_simple
'

sed "s/test_subset = .false./test_subset = .true./" input_base.nml > input.nml
test_expect_success "subset update" '
    mpirun -n 25 ../test_mpp_domains
'
sed "s/test_halosize_performance = .false./test_halosize_performance = .true./" input_base.nml > input.nml
test_expect_success "halosize performance" '
    mpirun -n 2 ../test_mpp_domains
'
sed "s/test_edge_update = .false./test_edge_update = .true./" input_base.nml > input.nml
test_expect_success "edge update" '
    mpirun -n 2 ../test_mpp_domains
'
sed "s/test_nonsym_edge = .false./test_nonsym_edge = .true./" input_base.nml > input.nml
test_expect_success "nonsym edge" '
    mpirun -n 2 ../test_mpp_domains
'
sed "s/test_performance = .false./test_performance = .true./" input_base.nml > input.nml
#If the system is Darwin or TRAVIS it will be skipped because it fails
test_expect_success "performance" '
    mpirun -n 6 ../test_mpp_domains
'
sed "s/test_global_sum = .false./test_global_sum = .true./" input_base.nml > input.nml
test_expect_success "global sum" '
    mpirun -n 2 ../test_mpp_domains
'
sed "s/test_cubic_grid_redistribute = .false./test_cubic_grid_redistribute = .true./" input_base.nml > input.nml
test_expect_success "cubic grid redistribute" '
    mpirun -n 6 ../test_mpp_domains
'
sed "s/test_boundary = .false./test_boundary = .true./" input_base.nml > input.nml
test_expect_success "boundary" '
    mpirun -n 6 ../test_mpp_domains
'
sed "s/test_adjoint = .false./test_adjoint = .true./" input_base.nml > input.nml
test_expect_success "adjoint" '
    mpirun -n 2 ../test_mpp_domains
'
sed "s/test_unstruct = .false./test_unstruct = .true./" input_base.nml > input.nml
test_expect_success "unstruct" '
    mpirun -n 2 ../test_mpp_domains
'
sed "s/test_group = .false./test_group = .true./" input_base.nml > input.nml
test_expect_success "group" '
    mpirun -n 2 ../test_mpp_domains
'
sed "s/test_interface = .false./test_interface = .true./" input_base.nml > input.nml
test_expect_success "interface" '
    mpirun -n 6 ../test_mpp_domains
'
sed "s/check_parallel = .false./check_parallel = .true./" input_base.nml > input.nml
test_expect_success "check_parallel" '
    mpirun -n 6 ../test_mpp_domains
'
sed "s/test_get_nbr = .false./test_get_nbr = .true./" input_base.nml > input.nml
test_expect_success "get nbr" '
    mpirun -n 8 ../test_mpp_domains
'

test_done
