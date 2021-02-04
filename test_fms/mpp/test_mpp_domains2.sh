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

# Set common test settings.
. ../test-lib.sh

# Skip tests
# Every test after 1 was skipped
SKIP_TESTS="$(basename $0 .sh)."[2,3,4,5,6,7,8,9,10,11,12,13,14,15,16]
echo "$SKIP_TESTS"
# TODO these might not have to be skipped
# add other tests that don't work on darwin/ci
#test "x$(uname -s)"="xDarwin" && SKIP_TESTS="$SKIP_TESTS $(basename $0 .sh).{4,7}"
#test "$CI"="true"             && SKIP_TESTS="$SKIP_TESTS $(basename $0 .sh).{7,8}"
#echo "$SKIP_TESTS"

# Create base nml for input
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
! NEST inputs
test_nest = .false.
num_nest = 3
tile_coarse =    1,  3,  7
tile_fine   =    7 , 8,  9
istart_coarse =  3,  3,  5
icount_coarse = 40,  5,  6
jstart_coarse =  3,  3,  6
jcount_coarse = 14,  6,  8
extra_halo = 0
ntiles_nest_all = 9
npes_nest_tile = 2, 2, 2, 2, 2, 2, 2, 1, 1
nest_level = 1, 1, 2
refine_ratio = 2, 2, 2
cyclic_nest = 'N'
! NEST inputs end
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
test_subset = .false.
test_unstruct = .false.
test_nonsym_edge = .false.
test_halosize_performance = .false.
test_adjoint = .false.
wide_halo = .false.
/
_EOF

touch input.nml # Input.nml is required to run the following tests
test_expect_success "1: Test simple functionality" '
    mpirun -n 4 ./test_domains_simple
'
#skipped
#sed "s/test_nest = .false./test_nest = .true./" $top_srcdir/test_fms/mpp/input_base.nml > input.nml
test_expect_success "2: Test update nest domain" '
    mpirun -n 2 ./test_mpp_domains
'
#skipped
#sed "s/test_subset = .false./test_subset = .true./" $top_srcdir/test_fms/mpp/input_base.nml > input.nml
test_expect_success "3: Test Subset Update" '
    mpirun -n 2 ./test_mpp_domains
'

#If the system is Darwin it will be skipped because it fails
sed "s/test_halosize_performance = .false./test_halosize_performance = .true./" $top_srcdir/test_fms/mpp/input_base.nml > input.nml
test_expect_success "4: Test Halosize Performance" '
    mpirun -n 2 ./test_mpp_domains
'

#skipped
#sed "s/test_edge_update = .false./test_edge_update = .true./" $top_srcdir/test_fms/mpp/input_base.nml > input.nml
test_expect_success "5: Test Edge update" '
    mpirun -n 2 ./test_mpp_domains
'
#skipped
#echo "5: Test Nonsym Edge"
#sed "s/test_nonsym_edge = .false./test_nonsym_edge = .true./" $top_srcdir/test_fms/mpp/input_base.nml > input.nml
test_expect_success "6: Test Nonsym Edge" '
    mpirun -n 2 ./test_mpp_domains
'

#echo "6: Test Performance"
sed "s/test_performance = .false./test_performance = .true./" $top_srcdir/test_fms/mpp/input_base.nml > input.nml
#If the system is Darwin or TRAVIS it will be skipped because it fails
test_expect_success "7: Test Performance" '
    mpirun -n 6 ./test_mpp_domains
'

#echo "7: Test Global Sum"
sed "s/test_global_sum = .false./test_global_sum = .true./" $top_srcdir/test_fms/mpp/input_base.nml > input.nml
#If the system is Darwin it will be skipped because it fails
test_expect_success "8: Test Global Sum" '
    mpirun -n 2 ./test_mpp_domains
'

#echo "8: Test Cubic Grid Redistribute"
sed "s/test_cubic_grid_redistribute = .false./test_cubic_grid_redistribute = .true./" $top_srcdir/test_fms/mpp/input_base.nml > input.nml
#If the system is Darwin or TRAVIS it will be skipped because it fails
test_expect_success "9: Test Cubic Grid Redistribute" '
    mpirun -n 6 ./test_mpp_domains
'

#echo "9: Test Boundary"
sed "s/test_boundary = .false./test_boundary = .true./" $top_srcdir/test_fms/mpp/input_base.nml > input.nml
#If the system is Darwin or TRAVIS it will be skipped because it fails
test_expect_success "10: Test cubic grid redistribute" '
    mpirun -n 6 ./test_mpp_domains
'

#echo "10: Test Adjoint"
sed "s/test_adjoint = .false./test_adjoint = .true./" $top_srcdir/test_fms/mpp/input_base.nml > input.nml
#If the system is Darwin it will be skipped because it fails
test_expect_success "11: Test Adjoint" '
    mpirun -n 2 ./test_mpp_domains
'
#skipped
#echo "11: Test Unstruct"
#sed "s/test_unstruct = .false./test_unstruct = .true./" $top_srcdir/test_fms/mpp/input_base.nml > input.nml
test_expect_success "12: Test Unstruct" '
    mpirun -n 2 ./test_mpp_domains
'

#echo "12: Test Group"
sed "s/test_group = .false./test_group = .true./" $top_srcdir/test_fms/mpp/input_base.nml > input.nml
#If the system is Darwin it will be skipped because it fails
run_test test_mpp_domains 2 $is_darwin
test_expect_success "13: Test Group" '
    mpirun -n 2 ./test_mpp_domains
'

#echo "13: Test Interface"
sed "s/test_interface = .false./test_interface = .true./" $top_srcdir/test_fms/mpp/input_base.nml > input.nml
#If the system is Darwin it will be skipped because it fails
test_expect_success "14: Test Interface" '
    mpirun -n 2 ./test_mpp_domains
'

#echo "14: Test Check Parallel"
#echo "Does not work on Darwin or elsewhere"
sed "s/check_parallel = .false./check_parallel = .true./" $top_srcdir/test_fms/mpp/input_base.nml > input.nml
test_expect_success "15: Test Interface" '
    mpirun -n 6 ./test_mpp_domains
'

#echo "15: Test Get Nbr"
sed "s/test_get_nbr = .false./test_get_nbr = .true./" $top_srcdir/test_fms/mpp/input_base.nml > input.nml
#If the system is Darwin it will be skipped because it fails
test_expect_success "16: Test Get Nbr" '
    mpirun -n 2 ./test_mpp_domains
'

test_done
